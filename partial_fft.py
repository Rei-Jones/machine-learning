def partial_fft(x):
    N = len(x)
    if N == 1: 
        return x, {}

    # Pad to next power of 2 if needed
    if (N & (N - 1)) != 0:
        pow_2 = 2**int(np.ceil(np.log2(N)))
        x = np.pad(x, (0, pow_2 - N), mode='constant')  
        N = len(x)

    # Bit reversal
    number = int(np.log2(N))
    bit_reversed = np.zeros(N, dtype=complex)
    for n in range(N):
        index = 0
        for a in range(number):
            if (n >> a) & 1:
                index |= (1 << (number - 1 - a))
        bit_reversed[index] = x[n]

    # Cooley-Tukey FFT
    for stage in range(1, number + 1): 
        size = 2**stage
        w_size = np.exp((-2j * np.pi) / size)
        for i in range(0, N, size):
            w = 1
          #Butterfly operation
            for j in range(size // 2):
                e = bit_reversed[i + j]
                o = w * bit_reversed[i + j + size // 2]
                bit_reversed[i + j] = e + o
                bit_reversed[i + j + size // 2] = e - o
                w *= w_size
    return bit_reversed, {}
