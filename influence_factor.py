import math
import numpy as np


def calculate_I1(B_prime, L_prime, H):
    M = L_prime / B_prime
    N = H / B_prime
    C1 = (1 + math.sqrt(M**2 + 1)) * math.sqrt(M**2 + N**2)
    C2 = M * (1 + math.sqrt(M**2 + N**2 + 1))
    C3 = (M + math.sqrt(M**2 + 1)) * math.sqrt(1 + N**2)
    C4 = M + math.sqrt(M**2 + N**2 + 1)
    I1 = (1 / math.pi) * (M * math.log(C1/C2) + math.log(C3/C4))

    return I1


def calculate_I2(B_prime, L_prime, H):
    M = L_prime / B_prime
    N = H / B_prime
    C1 = N  / ( 2 * math.pi)
    C2 = (np.arctan(M / (N * math.sqrt(M**2 + N**2 + 1))))
    I2 = C1 * C2

    return I2


def calculate_Is(B_prime, L_prime, H, nu):
    I1 = calculate_I1(B_prime, L_prime, H)
    I2 = calculate_I2(B_prime, L_prime, H)
    Is = I1 + ((1-2*nu)/(1-nu)) * I2

    return Is