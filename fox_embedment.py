import math
'''
Reference: Fox, E.N., 1948. The mean elastic settlement of a uniformly-loaded area at a depth
            below the ground surface, Proceedings of the 2nd International Conference, ISSMFE,
            Rotterdam, Vol. 1, pp. 129-132
'''


# Define a function to calculate the r factors
def calculate_r_factors(B, L, D):
    r = 2 * D
    r1 = math.sqrt(B ** 2 + r ** 2)
    r2 = math.sqrt(L ** 2 + r ** 2)
    r3 = math.sqrt(B ** 2 + L ** 2 + r ** 2)
    r4 = math.sqrt(B ** 2 + L ** 2)
    return r, r1, r2, r3, r4

# Define a function to calculate the beta factors
def calculate_beta_factors(nu):
    beta_1 = 3 - 4 * nu
    beta_2 = 5 - 12 * nu + 8 * nu ** 2
    beta_3 = -4 * nu * (1 - 2 * nu)
    beta_4 = -1 + 4 * nu - 8 * nu ** 2
    beta_5 = -4 * (1 - 2 * nu) ** 2
    return beta_1, beta_2, beta_3, beta_4, beta_5

# Define a function to calculate the Y factors
def calculate_Y_factors(B, L, r, r1, r2, r3, r4):
    Y_1 = B * math.log((r4 + L) / B) + L * math.log((r4 + B) / L) - (r4 ** 3 - B ** 3 - L ** 3) / (3 * B * L)
    Y_2 = B * math.log((r3 + L) / r1) + L * math.log((r3 + B) / r2) - (r3 ** 3 - r2 ** 3 - r1 ** 3 + r ** 3) / (3 * B * L)
    Y_3 = (r ** 2 / B) * math.log((L + r2) * r1 / ((L + r3) * r)) + (r ** 2 / L) * math.log((B + r1) * r2 / ((B + r3) * r))
    Y_4 = r ** 2 * (r1 + r2 - r3 - r) / (B * L)
    Y_5 = r * math.atan((B * L) / (r * r3))
    return Y_1, Y_2, Y_3, Y_4, Y_5

# Define a function to calculate the Fox Embedment Factor
def fox_factor(B, L, D, nu):
    r, r1, r2, r3, r4 = calculate_r_factors(B, L, D)
    beta_1, beta_2, beta_3, beta_4, beta_5 = calculate_beta_factors(nu)
    Y_1, Y_2, Y_3, Y_4, Y_5 = calculate_Y_factors(B, L, r, r1, r2, r3, r4)
    I_f = (beta_1*Y_1 + beta_2*Y_2 + beta_3*Y_3 + beta_4*Y_4 + beta_5*Y_5) / ((beta_1 + beta_2) * Y_1)
    return I_f


