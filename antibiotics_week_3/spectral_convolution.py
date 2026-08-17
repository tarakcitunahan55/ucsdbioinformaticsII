"""
SpectralConvolution
=================================="""

def read_experimental_spectrum(file):
    """Reads a single line of space separated integer masses."""
    with open(file, "r") as f:
        return list(map(int,f.readline().split())) #return [int(x) for x in f.readline().split()]


def spectral_convolution(spectrum):
    n=len(spectrum)
    convolution =[]
    for x in range(0,n-1):
        for y in range(x+1,n):
            diff=spectrum[y]-spectrum[x] #since the spectrum is given in increasing order all diff are >=0
            if diff!=0: #only count positive differences (not zero)
                convolution.append(diff)

    return convolution




if __name__ == "__main__":
    experimental_spectrum = read_experimental_spectrum("antibiotics_week_3/spectrum.txt")
    convolution=spectral_convolution(experimental_spectrum)
    print (*convolution)
