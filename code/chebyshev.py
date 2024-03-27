import matplotlib.pyplot as plt
import numpy as np
import math as math

pi_by_two = math.pi/2
forty_five_degrees = 0.7853981633974483 
methods = ['taylor', 'chebyshev-one-reduction', 'chebyshev-two-reductions']

def fatorial(number):
    if number == 0:
        return 1
    else:
        return number * fatorial(number-1)

def pow_degree(degree):
    return degree * degree

def angle_to_rad(degree):
    return degree * math.pi / 180

def taylor_even_series(degree):
    return 1 - pow(degree,2)/fatorial(2) + pow(degree,4)/fatorial(4) - pow(degree,6)/fatorial(6) + pow(degree,8)/fatorial(8) - pow(degree,10)/fatorial(10) + pow(degree,12)/fatorial(12)

def taylor_odd_series(degree):
    return degree - pow(degree,3)/fatorial(3) + pow(degree,5)/fatorial(5) - pow(degree,7)/fatorial(7) + pow(degree,9)/fatorial(9) - pow(degree,11)/fatorial(11)

def cos_taylor(degree):
    if degree<=forty_five_degrees and degree>=-forty_five_degrees:
        return taylor_even_series(degree)
    else:
        k = math.ceil((degree - forty_five_degrees)/pi_by_two)
        resulting_degree = degree - k*pi_by_two
        if k%4 == 0:
            return cos_taylor(resulting_degree)
        elif k%4 == 1:
            return -sin_taylor(resulting_degree)
        elif k%4 == 2:
            return -cos_taylor(resulting_degree)
        else:
            return sin_taylor(resulting_degree)

def sin_taylor(degree):
    if degree<=forty_five_degrees and degree>=-forty_five_degrees:
        return taylor_odd_series(degree)
    else:
        return cos_taylor(degree - pi_by_two)


def sin(degree): 
    if degree<=forty_five_degrees and degree>=-forty_five_degrees:
        y = pow_degree(degree)
        return degree*(3715891199/3715891200 + y * (-30965759/185794560 + y * (276479/33177600 + y * (-2879/14515200 + y * (13/4838400)))))
    else:
        return cos(degree - pi_by_two)


def sin2(degree): 
    if degree<=forty_five_degrees and degree>=-forty_five_degrees:
        y = pow_degree(degree)
        return degree*(116121589/116121600 + y * (-6193105/37158912 + y * (19343/2322432 + y * (-319/1658880))))
    else:
        return cos(degree - pi_by_two)



def cos(degree): 
    if degree<=forty_five_degrees and degree>=-forty_five_degrees:
        y = pow_degree(degree)
        return 980995276801/980995276800 + y * (-6812467201/13624934400 + y * (48660481/1167851520 + y * (-380161/273715200 + y * (503/20275200 + y * (-1/3548160)))))
    else:
        k = math.ceil((degree - forty_five_degrees)/pi_by_two)
        resulting_degree = degree - k*pi_by_two 
        if k%4 == 0:
            return cos(resulting_degree)
        elif k%4 == 1:
            return -sin(resulting_degree)
        elif k%4 == 2:
            return -cos(resulting_degree)
        else:
            return sin(resulting_degree)


def cos2(degree): 
    if degree<=forty_five_degrees and degree>=-forty_five_degrees:
        y = pow_degree(degree)
        return 12740198393/12740198400 + y * (-309657583/619315200 + y * (30965597/743178240 + y * (-138179/99532800 + y * (311/12902400))))
    else:
        k = math.ceil((degree - forty_five_degrees)/pi_by_two)
        resulting_degree = degree - k*pi_by_two 
        if k%4 == 0:
            return cos(resulting_degree)
        elif k%4 == 1:
            return -sin(resulting_degree)
        elif k%4 == 2:
            return -cos(resulting_degree)
        else:
            return sin(resulting_degree)

        
def plot_values(method):
    degree = []
    sines_errors = []
    cosines_errors = []

    actual_degree = -360

    if method == 'taylor':
        while actual_degree <= 360:
            degree.append(actual_degree)
            sines_errors.append(abs(math.sin(angle_to_rad(actual_degree)) - sin_taylor(angle_to_rad(actual_degree))))
            cosines_errors.append(abs(math.cos(angle_to_rad(actual_degree)) - cos_taylor(angle_to_rad(actual_degree))))

            actual_degree += 1
    elif method == 'chebyshev-one-reduction':
        while actual_degree <= 360:
            degree.append(actual_degree)
            sines_errors.append(abs(sin_taylor(angle_to_rad(actual_degree)) - sin(angle_to_rad(actual_degree))))
            cosines_errors.append(abs(cos_taylor(angle_to_rad(actual_degree)) - cos(angle_to_rad(actual_degree))))

            actual_degree += 1
    elif method == 'chebyshev-two-reductions':
        while actual_degree <= 360:
            degree.append(actual_degree)
            sines_errors.append(abs(sin_taylor(angle_to_rad(actual_degree)) - sin2(angle_to_rad(actual_degree))))
            cosines_errors.append(abs(cos_taylor(angle_to_rad(actual_degree)) - cos2(angle_to_rad(actual_degree))))

            actual_degree += 1
    
    return degree,sines_errors,cosines_errors

for method in methods:
    degree, sines_errors, cosines_errors = plot_values(method)

    plt.figure().set_size_inches(10, 10)

    plt.subplot(2, 1, 1)
    plt.plot(degree, sines_errors, color='blue')
    plt.title('Seno')
    plt.xlabel('Ângulo em graus')    
    plt.ylabel('Erro')
    plt.legend('Seno')
    plt.grid()

    plt.subplot(2, 1, 2)
    plt.plot(degree, cosines_errors, color='red')
    plt.title('Cosseno')
    plt.xlabel('Ângulo em graus')    
    plt.ylabel('Erro')
    plt.legend('Cosseno')
    plt.grid()

    plt.tight_layout()

    image = '../images/' + method + '.png'

    plt.savefig(image, dpi=300)