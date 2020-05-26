import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp
import math
import cmath
import os
import time

j = complex(0, 1)
wavelength = 663e-9
K = (2 * np.pi) / wavelength
delta_z = 0.01
arr_length = 51
arr_middle = int((arr_length + 1) / 2)
arr_half = int((arr_length - 1)/2)
gridsize = 0.001
pixel_pitch = gridsize / arr_length
length_fiber = 5
radius = 5e-6
max_number_of_modes = 10000                 # maximum l value used in characteristic equation root finding
n_core = 1.4700                             # 1.4450 for two mode output // 1.4650 for 20 mode output
n_clad = 1.4440
num_ap = np.sqrt(n_core**2 - n_clad**2)     # numerical aperture
v = K * radius * num_ap                     # waveguide parameter
u_max = K * radius * num_ap                 # maximum value of transverse wave number u, equal to the wave guide parameter v
u = np.linspace(0, u_max, 1001)             # array of u from 0 to u_max in steps of user input value
w = []
l = []


class Grids:

    def make_grid(self, x, y, fill):
        grid = np.full((x, y), fill)
        return grid

    def fill_grid_square(self, grid):
        rows = len(grid)
        columns = len(grid[0])
        grid[math.ceil(rows/2)-1-5:math.ceil(rows/2)-1+5, math.ceil(columns/2)-1-5:math.ceil(columns/2)-1+5] = 1
        #grid[0:101, 0:101] = 1
        return grid

    def fill_grid_cricle(self):
        vector = np.linspace(-50, 50, arr_length)
        radius = 10
        x, y = np.meshgrid(vector, vector)
        distance = np.sqrt(np.multiply(x, x) + np.multiply(y, y))
        a = np.where(distance <= radius)    # np.where will tell you where in the 2D array this condition is true
        b = np.where(distance > radius)
        distance[a] = 1
        distance[b] = 0
        return distance


# Calculation of the transverse wave number w
def w_calc(u_hold):
    w_hold = 0
    w_squared = v**2 - u_hold**2
    if w_squared < 0:
        w_hold = cmath.sqrt(w_squared)
    elif w_squared >= 0:
        w_hold = np.sqrt(w_squared)
    return w_hold


# characteristic equation with changing orders depending on the l value
def char_eq(x, y, order):
    solution = (x * sp.jv((order-1), x) / sp.jv(order, x)) + (y * sp.kv((order-1), y) / sp.kv(order, y))
    return solution


# Finds where the root is for array of characteristic equation points. This function only appends the roots that go from positive to negative
# If there are no roots left to be found the function will return 0, which will break the loop of root finding
# The function will only find up to one before the last value to avoid over running the array
def where_root_in_array(bessel):
    root_number_place = []
    sign = np.sign(bessel)
    for number in range(0, len(u)):
        if number + 1 == len(u):
            break
        elif sign[number] - sign[number + 1] == 2:
            root_number_place.append(number)
    if len(root_number_place) == 0:
        return 0
    else:
        return root_number_place


# Used in root finding for bisection method
def bisection(input1, input2):
    mid = (input1 + input2) / 2.0
    return mid


# Root finder to more accurately find a root between two values where we know a root exists
# Uses bisection method
def root_finder(u_hold1, u_hold2, order):
    tolerance = 5e-5
    iteration = 0
    itr_num = []
    bisect_mid = []
    bessel_diff = []
    maximum = 1000

    while iteration < maximum:
        w_hold1 = w_calc(u_hold1)
        result1 = char_eq(u_hold1, w_hold1, order)

        u_hold_mid = bisection(u_hold1, u_hold2)
        w_hold_mid = w_calc(u_hold_mid)
        result3 = char_eq(u_hold_mid, w_hold_mid, order)

        if (result3 * result1).real >= 0.0:
            u_hold1 = u_hold_mid
        else:
            u_hold2 = u_hold_mid
        true_root = bisection(u_hold1, u_hold2)

        if abs(true_root - u_hold_mid) < tolerance:
            # print(format(true_root, ".100f"))
            # plt.plot(itr_num, bessel_diff, "x-")
            # plt.xlabel("Iteration number")
            # plt.ylabel("Difference in characteristic equation output")
            # plt.show()
            #
            # plt.plot(itr_num, bisect_mid, "x-")
            # plt.xlabel("Iteration number")
            # plt.ylabel("Bisection mid value")
            # plt.show()
            return true_root

        diff = result3.real - result1.real
        bessel_diff.append(diff)
        bisect_mid.append(u_hold_mid)
        iteration += 1
        itr_num.append(iteration)
    return None


# Uses the two equations in finding fibre modes for core and cladding
# The output is a 2d grid electric field, where amplitude and phase can be found in the usual way
def fibre_mode_calc(u_hold, w_hold, l_hold):
    angle = np.zeros((arr_length, arr_length))
    dist_amp = np.zeros((arr_length, arr_length))
    f_grid = np.zeros((arr_length, arr_length), dtype=complex)

    vector = np.linspace(((-0.5*(arr_length-1))*pixel_pitch), ((0.5*(arr_length - 1))*pixel_pitch), arr_length)
    a = (0.25 * (arr_length - 1)) * pixel_pitch
    x, y = np.meshgrid(vector, vector)

    c_field = x + (j * y)
    for numrow in range(0, len(c_field)):
        for numcol in range(0, len(c_field)):
            angle[numrow][numcol] = np.arctan2(c_field[numrow][numcol].imag, c_field[numrow][numcol].real)
            dist_amp[numrow][numcol] = np.sqrt(c_field[numrow][numcol].real**2 + c_field[numrow][numcol].imag**2)
            if dist_amp[numrow][numcol] <= a:
                f_grid[numrow][numcol] = (sp.jv(l_hold, (u_hold * dist_amp[numrow][numcol]) / a) / sp.jv(l_hold, u_hold)) * np.exp(j * l_hold * angle[numrow][numcol])
            elif dist_amp[numrow][numcol] > a:
                f_grid[numrow][numcol] = (sp.kv(l_hold, (w_hold * dist_amp[numrow][numcol]) / a) / sp.kv(l_hold, w_hold)) * np.exp(j * l_hold * angle[numrow][numcol])
    return f_grid


# Calculation of the propagation constant
def beta_calc(u_root):
    beta_val = np.sqrt(K**2 * n_core**2 - (u_root / radius)**2)
    return beta_val


# Normalisation function for modes
def find_normalize_grid(grid, grid_cc):
    result_complex_num = 0
    for numrow in range(0, arr_length):
        for numcol in range(0, arr_length):
            result_complex_num += (grid[numrow][numcol] * grid_cc[numrow][numcol])
    normalization_factor = np.sqrt(1 / result_complex_num)
    return normalization_factor.real


# Checks if the mode has been normalised properly
def normalization_checker(grid):
    check = 0
    for numrow in range(0, arr_length):
        for numcol in range(0, arr_length):
            check += (grid[numrow][numcol].real**2 + grid[numrow][numcol].imag**2)  # sum of pixel intensities
    print(check)


# Function is used to check if the modes have been normalised properly, the output should always be 1.
def fiber_mode_product_checker(grid, conjugate_grid):
    check = 0
    for numrow in range(0, arr_length):
        for numcol in range(0, arr_length):
            check += grid[numrow][numcol] * conjugate_grid[numrow][numcol]
    print(check)


# Calculates the best version of the image that the fibre can propagate, essentially how well the image can be represented by the fibre
def image_reconstruction(complex_shift, mode_grid):
    reconstruct_grid = np.full((arr_length, arr_length), 0, dtype=complex)
    image_mode_product = complex_shift * mode_grid
    for numrow in range(0, arr_length):
        for numcol in range(0, arr_length):
            reconstruct_grid[numrow][numcol] += image_mode_product[numrow][numcol]
    return reconstruct_grid


# Calculates the scrambled image at the end of the fibre
def image_fiber_output(complex_num, mode_grid, beta_shift):
    reconstruct_grid = np.full((arr_length, arr_length), 0, dtype=complex)
    output_image = complex_num * mode_grid * beta_shift
    for numrow in range(0, arr_length):
        for numcol in range(0, arr_length):
            reconstruct_grid[numrow][numcol] += output_image[numrow][numcol]
    return reconstruct_grid


def angle_matrix():
    theta_x = []
    theta_y = []
    for x in range(0, arr_length):
        arg_x = (((arr_middle - 1) - x) / delta_z) * pixel_pitch
        theta_x.append(np.arcsin(arg_x))
    for y in range(0, arr_length):
        arg_y = (((arr_middle - 1) - y) / delta_z) * pixel_pitch
        theta_y.append(np.arcsin(arg_y))
    theta_X, theta_Y = np.meshgrid(theta_x, theta_y)
    K_x = np.multiply(K, np.arcsin(theta_X))
    K_y = np.multiply(K, np.arcsin(theta_Y))
    K_z = np.sqrt(np.multiply(K, K) - np.multiply(K_x, K_x) - np.multiply(K_y, K_y))
    return K_z  # theta_X


def amplitude(source):
    amplitude_output = np.sqrt(np.multiply(source.real, source.real) + np.multiply(source.imag, source.imag))
    return amplitude_output


def angular_spectrum_propagation(E_field, K_z, Z):
    e = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(E_field)))
    e_z = e * np.exp(j * K_z * Z)
    E_field_z = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(e_z)))
    return E_field_z


def angular_spectrum_propagation_3D(E_field, K_z, Z):
    e = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(E_field)))
    e_z = np.zeros((modes, arr_length, arr_length), dtype=complex)
    for x in range(0, modes, 1):
         e_z[x] = np.multiply(e[x], np.exp(j * K_z * Z))
    E_field_z = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(e_z)))
    return E_field_z


def modes_through_phase_plate(source, phase_plate):
    after_plate = np.zeros((modes, arr_length, arr_length), dtype=complex)
    for x in range(0, modes, 1):
        after_plate[x] = np.multiply(source[x], np.exp(np.multiply(j, phase_plate)))
    return after_plate


def normalise_array(array):
    a = np.multiply(array, array.conjugate())
    b = a.sum(axis=0)
    c = b.sum(axis=0)
    norm_array = np.divide(array, c)
    return norm_array


def conjugate_sum(mode_conjugates):
    sum = np.zeros((arr_length, arr_length), dtype=complex)
    for x in range(0, modes, 1):
        sum += mode_conjugates[x]
    return sum


def free_space_GS(source, target, K_z, pp_number):
    itteration = 0
    phase_plate = np.zeros((pp_number, arr_length, arr_length), dtype=complex)
    while itteration < 4:
        source_number = target_number = pp_number - 1
        m = 1
        n = 0
        A = target
        while n < pp_number:
            B = angular_spectrum_propagation_3D(source, K_z, (source_number * delta_z))
            C = np.multiply(B.conjugate(), A)
            J = conjugate_sum(C)
            phase_plate[source_number] = np.angle(J)
            D = modes_through_phase_plate(A, -phase_plate[source_number])
            source_number -= 1
            n += 1
            A = angular_spectrum_propagation_3D(D, K_z, -delta_z)
        E = modes_through_phase_plate(source, phase_plate[0])
        while m < pp_number:
            F = angular_spectrum_propagation_3D(E, K_z, delta_z)
            H = target
            if m < target_number:
                while target_number > m:
                    G = modes_through_phase_plate(H, -phase_plate[target_number])
                    H = angular_spectrum_propagation_3D(G, K_z, -delta_z)
                    target_number -= 1
            target_number = pp_number-1
            I = np.multiply(F.conjugate(), H)
            K = conjugate_sum(I)
            phase_plate[m] = np.angle(K)
            E = modes_through_phase_plate(F, phase_plate[m])
            m += 1
        itteration += 1
    print(itteration, "iterations run")
    A = modes_through_phase_plate(source, phase_plate[0])
    o = 1
    while o < pp_number:
        B = angular_spectrum_propagation_3D(A, K_z, delta_z)
        A = modes_through_phase_plate(B, phase_plate[o])
        o += 1
    approximate_field = A
    return approximate_field, phase_plate


def call_gauss_spot_2():
    spot = np.zeros((modes, arr_length, arr_length))
    final = np.zeros((arr_length, arr_length))
    vector = np.linspace(-arr_half, arr_half, arr_length)
    X, Y = np.meshgrid(vector, vector)
    sqrt_modes = np.sqrt(modes)
    row_number = round(sqrt_modes)
    majority_modes = round(modes / row_number)
    final_row = modes - ((row_number-1) * (majority_modes))
    if final_row > majority_modes:
        row_number += 1
        final_row -= majority_modes
    split_array_row_1 = round(arr_length / (majority_modes + 1))
    split_array_row_2 = round(arr_length / (final_row + 1))
    split_columns = round(arr_length / (row_number + 1))
    c1 = 1
    c2 = 1
    m = 0
    while c2 < row_number:
        while c1 <= majority_modes:
            spot[m] = two_d_gauss_spot(X, Y, (c1 * split_array_row_1) - arr_middle, (c2 * split_columns)-arr_middle)
            m += 1
            c1 += 1
        c1 = 1
        c2 += 1
    while c1 <= final_row:
        spot[m] = two_d_gauss_spot(X, Y, (c1 * split_array_row_2) - arr_middle, (c2 * split_columns)- arr_middle)
        m += 1
        c1 += 1
    n = 0
    while n < modes:
        final += spot[n]
        n+=1
    return spot


def two_d_gauss_spot(x, y, x_0, y_0):
    sigma_x = 1
    sigma_y = 1

    gauss_spot = np.exp(-((x - x_0)**2 / (2*sigma_x**2) + (y - y_0)**2 / (2 * sigma_y**2)))
    return gauss_spot


# *******************************
# *     APPLICATION START       *
# *******************************
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# *******************************
# *     FINDING FIBRE MODES     *
# *******************************
start_time = time.time()

modes_3d_original = []
modes_3d_scrambled = []
modes_counter = 0
G = Grids()
if os.path.exists("weighted_shift_complex_num.txt"):
    os.remove("weighted_shift_complex_num.txt")

img_file = open("face.txt", "r") # Input image is always 51x51, need ro find a way of excluding this
amplitude_grid_image = []
for lines in img_file:
    number_strings = lines.strip()
    numbers = [int(n) for n in number_strings]
    amplitude_grid_image.append(numbers)    # Amplitude of original image
phase_grid_image = G.make_grid(arr_length, arr_length, (np.pi/2))   # phase of original image.

reconstruct_img_complex_grid = np.full((arr_length, arr_length), 0, dtype=complex)
output_img_complex_grid = np.full((arr_length, arr_length), 0, dtype=complex)
product_test_mode = np.full((arr_length, arr_length), 0, dtype=complex)

# MODE FINDING
for num in range(0, len(u)):
    w_holder = w_calc(u[num])
    w.append(w_holder)      # array to hold all the transverse wave number w values

for num in range(0, max_number_of_modes):
    l.append(num)   # appending l values to an array from 0 up to a maximum number, this is used in the characteristic equation

for l_num in range(0, len(l)):
    roots = []  # array only holds u-roots
    bessel_out_real = char_eq(u, w, l[l_num]).real  # only need real values of the characteristic equation solution
    where_root = where_root_in_array(bessel_out_real)   # shows where the roots lie in the index of bessel output array
    if where_root == 0:
        print("No more roots to be found.") # breaks out of the loop as once there is no roots found, there will be no further roots found
        break
    else:
        for num in range(0, len(where_root)):
            index = where_root[num]     # using the array of the u roots to give the root finder two values to where the root will lie
            guess1 = u[index]
            guess2 = u[index + 1]
            roots.append(root_finder(guess1, guess2, l[l_num]))     # root finder finds the root to a more accurate value

    for count in range(0, len(roots)):
        w_val = w_calc(roots[count])    # corresponding w array for unique values of u
        beta = beta_calc(roots[count])  # propagation constant array

        print(l[l_num], count, roots[count], w_val, beta)   # l, m, u, v, beta

        mode_grid_fiber = fibre_mode_calc(roots[count], w_val, l[l_num])    # complex number 2d grid of the fibre mode

        mode_grid_fiber_cc = np.conjugate(mode_grid_fiber)
        normalization_factor = find_normalize_grid(mode_grid_fiber, mode_grid_fiber_cc)
        mode_grid_fiber_normalized = normalization_factor * mode_grid_fiber
        mode_grid_fiber_normalized_cc = np.conjugate(mode_grid_fiber_normalized)

        modes_3d_original.append(mode_grid_fiber_normalized)

        amp_grid_fiber = np.sqrt(mode_grid_fiber_normalized.real**2 + mode_grid_fiber_normalized_cc.imag**2)
        phase_grid_fiber = np.arctan2(mode_grid_fiber_normalized.imag, mode_grid_fiber_normalized.real)

        # SHOW THE AMPLITUDE AND MODES AS THEY ARE FOUND
        # plt.imshow(amp_grid_fiber, cmap='hot', interpolation='nearest')                   # Plots the amplitude and phase of each mode as they are found
        # plt.title("Amplitude graph")
        # cax = plt.axes([0.85, 0.1, 0.075, 0.8])
        # plt.colorbar(cax=cax)
        # plt.show()
        #
        # plt.imshow(phase_grid_fiber, cmap='hot', interpolation='nearest')
        # plt.title("Phase graph")
        # cax = plt.axes([0.85, 0.1, 0.075, 0.8])
        # plt.colorbar(cax=cax)
        # plt.show()

        #   IMAGE PROCESSING
        mode_normalised_cc_amp = np.sqrt(mode_grid_fiber_normalized_cc.real ** 2 + mode_grid_fiber_normalized_cc.imag ** 2)     # cc = complex conjugate
        mode_normalised_cc_phase = np.full((arr_length, arr_length), 0, dtype=float)
        for numrow in range(0, arr_length):
            for numcol in range(0, arr_length):
                mode_normalised_cc_phase[numrow][numcol] = math.atan2(mode_grid_fiber_normalized_cc[numrow][numcol].imag, mode_grid_fiber_normalized_cc[numrow][numcol].real)

        mode_grid_cc_polar = mode_normalised_cc_amp * np.exp(j * mode_normalised_cc_phase)
        image_grid_polar = amplitude_grid_image * np.exp(j * phase_grid_image)
        mode_image_total = mode_grid_cc_polar * image_grid_polar

        mode_image_sum_complex_num = 0
        for numrow in range(0, arr_length):
            for numcol in range(0, arr_length):
                mode_image_sum_complex_num += mode_image_total[numrow][numcol]
        with open("weighted_shift_complex_num.txt", "a+") as complex_num_shift_file:
            complex_num_shift_file.write(str(mode_image_sum_complex_num) + "\n")

        # IMAGE RECONSTRUCTION
        reconstruct_img_complex_grid += image_reconstruction(mode_image_sum_complex_num, mode_grid_fiber_normalized)    # adding all modes together to reform the image

        # IMAGE FIBER OUTPUT
        beta_phase_shift = np.exp(j * length_fiber * beta)
        output_img_complex_grid += image_fiber_output(mode_image_sum_complex_num, mode_grid_fiber_normalized, beta_phase_shift)     # adding all scrambled modes together

        modes_3d_scrambled.append(mode_grid_fiber_normalized * beta_phase_shift)
        modes_counter += 1

modes = modes_counter

reconstruct_img_amp = np.sqrt(reconstruct_img_complex_grid.real**2 + reconstruct_img_complex_grid.imag**2)
reconstruct_img_phase = np.full((arr_length, arr_length), 0, dtype=float)
for numrow in range(0, arr_length):
    for numcol in range(0, arr_length):
        reconstruct_img_phase[numrow][numcol] = math.atan2(reconstruct_img_complex_grid[numrow][numcol].imag, reconstruct_img_complex_grid[numrow][numcol].real)

output_img_amp = np.sqrt(output_img_complex_grid.real**2 + output_img_complex_grid.imag**2)     # output image amplitude
output_img_phase = np.full((arr_length, arr_length), 0, dtype=float)    # output image phase
for numrow in range(0, arr_length):
    for numcol in range(0, arr_length):
        output_img_phase[numrow][numcol] = math.atan2(output_img_complex_grid[numrow][numcol].imag, output_img_complex_grid[numrow][numcol].real)

# Graphs to show the fibre output and reconstructions
# plt.imshow(reconstruct_img_amp, cmap='hot', interpolation='nearest')
# cax = plt.axes([0.85, 0.1, 0.075, 0.8])
# plt.colorbar(cax=cax)
# plt.show()
# plt.imshow(reconstruct_img_phase, cmap='hot', interpolation='nearest')
# cax = plt.axes([0.85, 0.1, 0.075, 0.8])
# plt.colorbar(cax=cax)
# plt.show()
# plt.imshow(output_img_amp, cmap='hot', interpolation='nearest')
# cax = plt.axes([0.85, 0.1, 0.075, 0.8])
# plt.colorbar(cax=cax)
# plt.show()
# plt.imshow(output_img_phase, cmap='hot', interpolation='nearest')
# cax = plt.axes([0.85, 0.1, 0.075, 0.8])
# plt.colorbar(cax=cax)
# plt.show()


# *******************************
# *             MPLC            *
# *******************************

source_field = np.zeros((modes, arr_length, arr_length), dtype=complex)     # Here we are allocating arrays with the correct dimensions and data type
target_field = np.zeros((modes, arr_length, arr_length), dtype=complex)     # for use later on in the programme.
origonal_field = np.zeros((modes, arr_length, arr_length), dtype=complex)

for w in range(0, modes):
    origonal_field[w] = modes_3d_original[w]
    source_field[w] = modes_3d_scrambled[w]

target_field = call_gauss_spot_2()

# for y in range(0, modes):
#     plt.imshow(amplitude(source_field[y]), cmap='hot', interpolation='nearest')   # Plots of the amplitude and phase of the modes after they have propagated through fibre.
#     plt.title("Mode %s Amplitude After Fibre Propagation" % y)
#     plt.xlabel('Length %sm' % gridsize)
#     plt.ylabel('Length %sm' % gridsize)
#     plt.show()
#     plt.imshow(np.angle(source_field[y]), cmap='hot', interpolation='nearest', vmin=-np.pi, vmax=np.pi)
#     plt.title("Mode %s Phase After Fibre Propagation" % y)
#     plt.xlabel('Length %sm' % gridsize)
#     plt.ylabel('Length %sm' % gridsize)
#     plt.show()

angular_K_z = angle_matrix()
output_of_first_plates = np.zeros((modes, arr_length, arr_length), dtype=complex)   # Here we are allocating arrays with the correct dimaensions and data type
output_of_second_plates = np.zeros((modes, arr_length, arr_length), dtype=complex)  # for use later on in the programme.
output_of_third_plates = np.zeros((modes, arr_length, arr_length), dtype=complex)
fixed_phase_Gaussians = np.zeros((modes, arr_length, arr_length), dtype=complex)

output_of_first_plates, phase_plates_1 = free_space_GS(source_field, target_field, angular_K_z, 100)

gaussian_sum = np.zeros((arr_length, arr_length), dtype=complex)

for w in range(0, modes):
    gaussian_sum += output_of_first_plates[w]

plt.imshow(amplitude(gaussian_sum), cmap='hot', interpolation='nearest')    # Outputs all of the spatially seperated Gaussian spots onto one plot.
plt.title("Guassian Spots Amplitudes")
plt.xlabel('Length %sm' %gridsize)
plt.ylabel('Length %sm' %gridsize)
plt.show()
plt.imshow(np.angle(gaussian_sum), cmap='hot', interpolation='nearest', vmin=-np.pi, vmax=np.pi)
plt.title("Guassian Spots Phase")
plt.xlabel('Length %sm' %gridsize)
plt.ylabel('Length %sm' %gridsize)
plt.show()


# By the "output_of_third_plates" line we have used the 'free_space_GS' three times to create 3 different sets of phase plates that make up one system.
# The simulation also obtains the final output of each restored mode which is plotted in the next few lines.
# The first set of plates is responsible for arranging the modes into Gaussian spots.
# The second set, which is a single plate, is responsible from restoring the phase of each mode to its original form before it propagated through the fibre.
# The Final set of plates re-sorts the Guassian spots into the original image.

fixed_phase_Gaussians = np.multiply(amplitude(output_of_first_plates), np.exp(np.multiply(j, np.angle(origonal_field))))
output_of_second_plates, phase_plates_2 = free_space_GS(output_of_first_plates, fixed_phase_Gaussians, angular_K_z, 1)
output_of_third_plates, phase_plates_3 = free_space_GS(output_of_second_plates, origonal_field, angular_K_z, 100)

# for m in range(0, modes):
#     plt.imshow(amplitude(output_of_third_plates[m]), cmap='hot', interpolation='nearest')
#     plt.title("Restored Mode %s Amplitude" % m)
#     plt.xlabel('Length %sm' % gridsize)   # The plotted output here shows the apmlitude and phase of each mode after
#     plt.ylabel('Length %sm' % gridsize)   # it has been restored in the MPLC system.
#     plt.show()
#     plt.imshow(np.angle(output_of_third_plates[m]), cmap='hot', interpolation='nearest', vmin=-np.pi, vmax=np.pi)
#     plt.title("Restored Mode %s Phase" % m)
#     plt.xlabel('Length %sm' % gridsize)
#     plt.ylabel('Length %sm' % gridsize)
#     plt.show()

# IMAGE UNSCRAMBLING USING CREATED PHASE PLATES

# plt.imshow(amplitude_grid_image, cmap='hot', interpolation='nearest') # This section plots the amplitude and phase of the earlier defined field that
# plt.title("Original Smiley Face Amplitude")                           # we are propagating through our fibre. The First two plots will show the amplitude
# plt.xlabel('Length %sm' % gridsize)                                   # and phase of the field before propagating, the last two will show the apmplitude
# plt.ylabel('Length %sm' % gridsize)                                   # and phase after it has travelled through the fibre and become scrambled.
# plt.show()
# plt.imshow(phase_grid_image, cmap='hot', interpolation='nearest', vmin=-np.pi, vmax=np.pi)
# plt.title("Original Smiley Face Phase")
# plt.xlabel('Length %sm' % gridsize)
# plt.ylabel('Length %sm' % gridsize)
# plt.show()
# plt.imshow(output_img_amp, cmap='hot', interpolation='nearest')
# plt.title("Scrambled Smiley Face Amplitude")
# plt.xlabel('Length %sm' % gridsize)
# plt.ylabel('Length %sm' % gridsize)
# plt.show()
# plt.imshow(output_img_phase, cmap='hot', interpolation='nearest', vmin=-np.pi, vmax=np.pi)
# plt.title("Scrambled Smiley Face Phase")
# plt.xlabel('Length %sm' % gridsize)
# plt.ylabel('Length %sm' % gridsize)
# plt.show()

A = np.multiply(output_img_amp, np.exp(np.multiply(j, output_img_phase)))   # the code here propagates the scrambled field through each phase plate that
B = np.multiply(A, np.exp(np.multiply(j, phase_plates_1[0])))               # was created in the earlier program. The amplitude and phase is plotted at key
for x in range(1, 100):                                                     # intervals as the field passes through the MPLC. These key points are; after
    C = angular_spectrum_propagation(B, angular_K_z, delta_z)               # the input image has had its modes spacially divided into gaussian spots, and
    B = np.multiply(C, np.exp(np.multiply(j, phase_plates_1[x])))           # after the spots have been restored to the form of the orignal image.
D = np.multiply(B, np.exp(np.multiply(j, phase_plates_2[0])))

plt.imshow(amplitude(D), cmap='hot', interpolation='nearest')
plt.title("Unscrambled Smiley Face Amplitude")
plt.xlabel('Length %sm' % gridsize)
plt.ylabel('Length %sm' % gridsize)
plt.show()
plt.imshow(np.angle(D), cmap='hot', interpolation='nearest', vmin=-np.pi, vmax=np.pi)
plt.title("Unscrambled Smiley Face Phase")
plt.xlabel('Length %sm' % gridsize)
plt.ylabel('Length %sm' % gridsize)
plt.show()
E = np.multiply(D, np.exp(np.multiply(j, phase_plates_3[0])))
for y in range(1, 100):
    F = angular_spectrum_propagation(E, angular_K_z, delta_z)
    E = np.multiply(F, np.exp(np.multiply(j, phase_plates_3[y])))

plt.imshow(amplitude(E), cmap='hot', interpolation='nearest')
plt.title("Unscrambled Smiley Face Amplitude")
plt.xlabel('Length %sm' % gridsize)
plt.ylabel('Length %sm' % gridsize)
# cax = plt.axes([0.85, 0.1, 0.075, 0.8])
# plt.colorbar(cax=cax)
plt.show()
plt.imshow(np.angle(E), cmap='hot', interpolation='nearest', vmin=-np.pi, vmax=np.pi)
plt.title("Unscrambled Smiley Face Phase")
plt.xlabel('Length %sm' % gridsize)
plt.ylabel('Length %sm' % gridsize)
# cax = plt.axes([0.85, 0.1, 0.075, 0.8])
# plt.colorbar(cax=cax)
plt.show()

print("--- %s seconds to execute ---" % (time.time() - start_time))     # Program finished by outputting how much time was required to run the entire simulation.

