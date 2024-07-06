def get_excitation_degree(I1, I2):
    degree = 0
    for i in range(len(I1)):
        degree += bin(I1[i] ^ I2[i]).count('1')
    return int(degree / 2)