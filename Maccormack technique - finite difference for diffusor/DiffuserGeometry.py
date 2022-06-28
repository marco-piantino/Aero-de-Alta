import numpy as np

def totalLength(L_conv: float, L_t: float, L_div: float):
    L = L_conv + L_t + L_div
    return L


def radius(L_conv: float, L_t: float, L_div: float, r_0: float, r_t: float, r_f: float, x: float):
    radius = 0
    L = totalLength(L_conv, L_t, L_div)
    if x < L_conv:
        radius = r_0 + (r_t - r_0)*x/L_conv
    elif x < L_conv + L_t:
        radius = r_t
    elif x <= L:
        radius = r_t + (r_f - r_t)*(x - (L_conv + L_t))/L_div
    return radius
