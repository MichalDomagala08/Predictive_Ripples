

from email.mime import image
from xml.parsers.expat import model

from scipy.spatial.distance import cosine

import torch
import torch.nn.functional as F
from PIL import Image
from torchvision import models
from torchvision.models import resnet50, ResNet50_Weights
import pandas as pd
from pathlib import Path
from PIL import Image
import numpy as np
import cv2
import numpy as np

def prepare_model(model,model_name = 'resnet50'):

    model.eval()  # switch to inference mode (e.g., deactivate dropout and batch normalization layers)

    # freeze all parameters to avoid computing gradients (alhough we are not training)
    for p in model.parameters():
        p.requires_grad = False

    # Hooks - function that will be called during the forward pass to capture the output of specific layers
    # It saves activation inside the dictionairy that is empty. It saves module( layer3, layer4) output to the dictionary with the name of the layer as key.
    # inp - what was inputed and out - what was outputed (i.e. our activations!) - we take only ouptut ofc 
    feats = {}
    def make_hook(name):
        def hook(module, inp, out):
            feats[name] = out
        return hook

    if model_name == 'resnet50':
        model.layer3.register_forward_hook(make_hook("layer3")) # Register  hook for said layers
        model.layer4.register_forward_hook(make_hook("layer4"))
        sizes = [model.layer4[-1].conv3.out_channels, model.layer3[-1].conv3.out_channels]

    elif model_name == 'dino':
        model.blocks[-1].register_forward_hook(make_hook("block_last"))   # semantyka
        model.blocks[-4].register_forward_hook(make_hook("block_mid"))    # przestrzenno-tek
        sizes = [model.blocks[-1].norm2.normalized_shape[0], model.blocks[-1].norm2.normalized_shape[0]]

    return model,feats,sizes

def get_image_repr(image, model, tmf,feats,model_name = 'resnet50'):
    # Preprocessing:
    # - input either weights or tfm object (for dino)

    # Apply preprocessing used for this particular weight configuration to our image
    device = next(model.parameters()).device
    x = tmf(Image.fromarray(image)).unsqueeze(0).to(device, non_blocking=True)

        
    # Inference mode! Turn off the gradient coputation for memory efficiency
    with torch.inference_mode():
        model(x)

    # Get representations and global pool and normalize them:
    if model_name == 'resnet50':
        rep_contr = F.normalize(feats["layer3"].mean(dim=(2, 3)), dim=1).squeeze(0) 
        rep_main  = F.normalize(feats["layer4"].mean(dim=(2, 3)), dim=1).squeeze(0)  

    elif model_name == 'dino':
        B   = x.shape[0]
        ps  = model.patch_size                       # 14
        h, w = x.shape[2] // ps, x.shape[3] // ps    # 518/14 = 37 -> siatka 37x37
 
        # tokeny (B, 1+N, C) -> odciecie CLS -> mapa (B, C, h, w)
        fmap_last = feats["block_last"][:, 1:, :].reshape(B, h, w, -1).permute(0, 3, 1, 2)
        fmap_mid  = feats["block_mid"] [:, 1:, :].reshape(B, h, w, -1).permute(0, 3, 1, 2)

        rep_main  = F.normalize(fmap_last.mean(dim=(2, 3)), dim=1).squeeze(0)
        rep_contr = F.normalize(fmap_mid.mean(dim=(2, 3)),  dim=1).squeeze(0)

    return rep_main, rep_contr

def pixels_within_radius(image_4d, image_index, x, y, radius):
    image = image_4d[image_index]
    height, width = image.shape[:2]
    x, y, radius = int(x), int(y), int(radius)
    y_grid, x_grid = np.ogrid[:height, :width]
    inside = (x_grid - x) ** 2 + (y_grid - y) ** 2 <= radius ** 2
    pixel_coordinates = np.column_stack(np.where(inside))[:, ::-1]
    pixel_values = image[inside]
    return pixel_coordinates, pixel_values


def reconstruct_image(image_4d, pixel_coordinates, pixel_values, image_index):
    reconstructed = np.zeros_like(image_4d[image_index])
    x, y = np.astype(pixel_coordinates[:, 0],'int'), np.astype(pixel_coordinates[:, 1],'int')
    reconstructed[y, x] = pixel_values
    return reconstructed
       
from scipy.ndimage import distance_transform_edt


def make_alpha_mask(shape, pixel_coordinates, feather=8):
    """
    Maska alfa dla SUMY wszystkich dotychczasowych fiksacji.

    alpha = 1      -> srodek dyskow: oryginalny obraz bez zmian
    1 > alpha > 0  -> pierscien `feather` px na zewnatrz krawedzi
    alpha = 0      -> tlo

    Nachodzace sie dyski lacza sie w jeden region - feather biegnie
    po zewnetrznym konturze sumy, nie po kazdym dysku osobno.
    """
    mask = np.zeros(shape, dtype=bool)
    x = pixel_coordinates[:, 0].astype(int)
    y = pixel_coordinates[:, 1].astype(int)
    mask[y, x] = True

    # odleglosc kazdego piksela tla od najblizszego odsłoniętego piksela
    dist = distance_transform_edt(~mask)

    # w srodku dist=0 -> alpha=1; liniowy spadek do 0 przy dist=feather
    return np.clip(1.0 - dist / feather, 0.0, 1.0)


def reconstruct_image_2(image_4d, pixel_coordinates, pixel_values, image_index,
                      feather=8, background=None):
    """
    Sklada obraz wejsciowy dla modelu przez maske alfa.
    `pixel_values` nie jest juz potrzeba (piksele bierzemy z oryginalu tam,
    gdzie alpha=1) - zostaje w sygnaturze, zeby nie zmieniac wywolania.
    """
    image = image_4d[image_index].astype(np.float32)

    if len(pixel_coordinates) == 0:
        alpha = np.zeros(image.shape[:2], dtype=np.float32)
    else:
        alpha = make_alpha_mask(image.shape[:2], pixel_coordinates, feather)

    # Neutralne tlo = srednia CALEGO obrazu. Celowo stale w obrebie proby:
    # tlo musi byc identyczne na kazdym etapie, inaczej zmiany podobieństwa
    # moglyby byc napędzane przez dryf tła, a nie przez nową treść.
    if background is None:
        background = image.reshape(-1, 3).mean(axis=0)
    background = np.asarray(background, dtype=np.float32)

    out = alpha[..., None] * image + (1.0 - alpha[..., None]) * background
    return np.clip(out, 0, 255).astype(np.uint8)