epsilon = 1e-6

weights_0 = 0.1
weights_1 = 0.6
weights_2 = 0.3


def weno5_nonlinear_weights(value):
    v0 = value[0]
    v1 = value[1]
    v2 = value[2]
    v3 = value[3]
    v4 = value[4]

    beta00 = v0 - 2 * v1 + v2
    beta01 = v0 - 4 * v1 + 3 * v2
    beta0 = 13 / 12 * beta00**2 + 1 / 4 * beta01**2

    beta10 = v1 - 2 * v2 + v3
    beta11 = v1 - v3
    beta1 = 13 / 12 * beta10**2 + 1 / 4 * beta11**2

    beta20 = v2 - 2 * v3 +  v4
    beta21 = 3 * v2 - 4 * v3 + v4
    beta2 = 13 / 12 * beta20**2 + 1 / 4 * beta21**2

    alpha0 = weights_0 / (beta0 + epsilon)**2 
    alpha1 = weights_1 / (beta1 + epsilon)**2
    alpha2 = weights_2 / (beta2 + epsilon)**2

    w0 = alpha0 / (alpha0 + alpha1 + alpha2)
    w1 = alpha1 / (alpha0 + alpha1 + alpha2)
    w2 = alpha2 / (alpha0 + alpha1 + alpha2)
    
    return w0, w1, w2

def weno5_reconstruction(value):
    v0 = value[0]
    v1 = value[1]
    v2 = value[2]
    v3 = value[3]
    v4 = value[4]
    
    w0, w1, w2 = weno5_nonlinear_weights(value)
    
    term0 = w0 * (1/3 * v0 - 7/6 * v1 + 11/6 * v2)
    term1 = w1 * (-1/6 * v1 + 5/6 * v2 + 1/3 * v3)
    term2 = w2 * (1/3 * v2 + 5/6 * v3 - 1/6 *v4)
    
    return term0 + term1 + term2
    
