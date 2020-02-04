import numpy as np
from numpy.random import *
import matplotlib.pyplot as plt
from numpy.fft import *

""" Ce script genere des interferogrammes tels qu'obtenus avec un interferometre
de Michelson dans le but d'etudier la transformée de Fourier et de comprendre 
comment la resolution spectrale est déterminée.
"""

def read_from_file(filename):
	file = np.genfromtxt(filename, skip_header=18, skip_footer=1)
	x = file[:, 1]
	y = file[:, 2]
	return (x, y)  # J: J'ai chamgé la fonction np.loadtxt par la fonction np.genfromtxt (plus efficace)

def plot_ocean_optics(filename, title, left, right):
	file = np.genfromtxt(filename)
	plt.figure(figsize=(10,7))
	#plt.title("{}".format(title))
	plt.xlabel("Longueur d'onde (nm)")
	plt.plot(file[:,0],file[:,1])
	plt.xlim(left=left, right=right)

plot_ocean_optics("autre_lampe_sodium_faible.txt", "Ocean Optics lampe sodium faible", 400, 700)
plot_ocean_optics("autre_lampe_sodium_forte.txt", "Ocean Optics lampe sodium forte", 400, 900)
plot_ocean_optics("autre_Lazer_He_Ne.txt", "Ocean Optics laser He-Ne", 400, 700)
plot_ocean_optics("autre_lumiere_blanche.txt", "Ocean Optics lumiere blanche", 200, 900)

def generateHeNeInterferogram(xMin, xMax, N):
	""" Genere un tableau de N valeurs equidistantes enntre xMin et xMax.
	Ensuite, genere un tableau de N valeurs qui representent un interferogramme
	d'un laser He-Ne a 0.6328 microns. On ajoute du bruit pour rendre le tout
	plus realiste.
	"""
	dx = (xMax - xMin)/N
	x = np.linspace(xMin, xMax, N)
	noise = random(len(x))*0.05
	y = 1+np.cos(2 * np.pi / 0.6328 * x)+noise
	return (x,y)

def generateWhiteLightInterferogram(xMin, xMax, N):
	""" Genere un tableau de N valeurs equidistantes enntre xMin et xMax.
	Ensuite, genere un tableau de N valeurs qui representent un interferogramme
	d'une source blanche visible. On ajoute du bruit pour rendre le tout
	plus realiste.
	"""
	dx = (xMax - xMin)/N
	x = np.linspace(xMin, xMax, N)
	noise = random(len(x))*0.05
	k1 = 1/0.4
	k2 = 1/0.8
	y = 1+np.exp(-x*x/4)*(np.sin(2 * np.pi * (k1+k2)*x/2)/x * np.sin(2 * np.pi * (k1-k2)*x/2)+ noise)
	return (x,y)

def fourierTransformInterferogram(x,y):
	""" A partir du tableau de valeurs Y correspondant a l'abscisse X,
	la transformée de Fourier est calculée et l'axes des fréquences (f en
	µm^-1) et des wavelengths (1/f en microns) est retournée.

	Le spectre est un ensemble de valeurs complexes pour lesquelles l'amplitude
	et la phase sont pertinentes: l'ordre des valeurs commence par la valeur DC (0)
	et monte jusqu'a f_max=1/2/∆x par resolution de ∆f = 1/N/∆x. A partir de la
	(N/2) ieme valeur, la frequence est negative jusqu'a -∆f dans la N-1 case.
	Voir
	https://github.com/dccote/Enseignement/blob/master/HOWTO/HOWTO-Transformes%20de%20Fourier%20discretes.pdf
	"""
	spectrum = fft(y)
	dx = 2*(x[1]-x[0]) # on obtient dx, on suppose equidistant
	N = len(x)     # on obtient N directement des données
	frequencies = fftfreq(N, dx) # Cette fonction est fournie par numpy
	wavelengths = 1/frequencies-float(200e-9)  # Les fréquences en µm^-1 sont moins utiles que lambda en µm
	return (wavelengths, frequencies, spectrum)



def plotCombinedFigures(x, y, w, s, title="", left=400, right=800):
	""""
	On met l'interferogramme et le spectre sur la meme page.
	"""
	fig, (axes, axesFFT) = plt.subplots(2,1,figsize=(10, 7))
	axes.plot(x, y, '-m')
	axes.set_title("Interferogramme")
	axesFFT.plot(w*1000, abs(s), '.-m')
	axesFFT.set_xlim(left=left, right=right)
	axesFFT.set_xlabel("Longueur d'onde [nm]")
	axesFFT.set_title(title)


# Basse resolution
#(x,y) = generateHeNeInterferogram(xMin=-10, xMax=10, N=200) # en microns
#(w, f, s)  = fourierTransformInterferogram(x,y)
#df = f[1]-f[0]
#dl = 0.6328*0.6328*df*1000 # x 1000 pour nm
#plotCombinedFigures(x,y,w,s,left=632.8-5*dl, right=632.8+5*dl, title="Spectre He-Ne basse resolution {0:0.2f} nm".format(dl))

# Haute resolution
# Resolution ∆f = 1/(200 µm * 2000)
# Resolution @ 632.8 nm : ∆lambda = 632.8^2 * ∆f car ∆lambda/lambda = ∆f/f.
#(x,y) = generateHeNeInterferogram(xMin=-100, xMax=100, N=2000) # en microns
#(w, f, s) = fourierTransformInterferogram(x,y)
#df = f[1]-f[0]
#dl = 0.6328*0.6328*df*1000
#plotCombinedFigures(x,y,w,s,left=632.8-5*dl, right=632.8+5*dl, title="Spectre He-Ne haute resolution {0:0.2f} nm".format(dl))


# Tres haute resolution
# Resolution ∆f = 1/(2000 µm * 2000)
# Resolution @ 632.8 nm : ∆lambda = 632.8^2 * ∆f
#(x,y) = generateHeNeInterferogram(xMin=-1000, xMax=1000, N=20000) # en microns
#(w, f, s) = fourierTransformInterferogram(x,y)
#df = f[1]-f[0]
#dl = 0.6328*0.6328*df*1000
#plotCombinedFigures(x,y,w,s,left=632.8-5*dl, right=632.8+5*dl, title="Spectre He-Ne tres haute resolution {0:0.2f} nm".format(dl))


# Hyper haute resolution
# Resolution ∆f = 1/(20000 µm * 20000)
# Resolution @ 632.8 nm : ∆lambda = 632.8^2 * ∆f
#(x,y) = generateHeNeInterferogram(xMin=-10000, xMax=10000, N=200000) # en microns
#(w, f, s) = fourierTransformInterferogram(x,y)
#df = f[1]-f[0]
#dl = 0.6328*0.6328*df*1000
#plotCombinedFigures(x,y,w,s,left=632.8-5*dl, right=632.8+5*dl, title="Spectre He-Ne hyper haute resolution {0:0.2f} nm".format(dl))


# Spectre de lumiere blanche
# Resolution ∆f = 1/(20000 µm * 20000)
# Resolution @ 500 nm : ∆lambda = 500^2 * ∆f
#(x,y) = generateWhiteLightInterferogram(xMin=-100, xMax=100, N=20000) # en microns
#(w, f, s)  = fourierTransformInterferogram(x,y)
#df = f[1]-f[0]
#dl = 0.500*0.500*df*1000 # resolution autour de 0.500 µm en nm
#plotCombinedFigures(x,y,w,s,left=0, right=2000, title="Spectre lumiere blanche, resolution {0:0.2f} nm".format(dl))

(x, y) = read_from_file("laser_he_ne_0.txt")
(w, f, s) = fourierTransformInterferogram(x,y)
df = f[1] - f[0]
dl = 0.500*0.500*df*1000
#plotCombinedFigures(x, y, w, s, left=0, right=1000, title="Laser HeNe, resolution {0:0.2f} nm".format(dl))

(x, y) = read_from_file("laser_he_ne_1.txt")
(w, f, s) = fourierTransformInterferogram(x,y)
df = f[1] - f[0]
dl = 0.500*0.500*df*1000
#plotCombinedFigures(x, y, w, s, left=0, right=800, title="Laser HeNe, resolution {0:0.2f} nm".format(dl))

(x, y) = read_from_file("laser_he_ne_3.txt")
(w, f, s) = fourierTransformInterferogram(x,y)
df = f[1] - f[0]
dl = 0.500*0.500*df*1000
plotCombinedFigures(x, y, w, s, left=300, right=800, title="Laser He-Ne, resolution {0:0.2f} nm".format(dl))

(x, y) = read_from_file("lum_blanche_0.txt")
(w, f, s) = fourierTransformInterferogram(x,y)
df = f[1] - f[0]
dl = 0.500*0.500*df*1000
#plotCombinedFigures(x, y, w, s, left=0, right=800, title="Lumière blanche, resolution {0:0.2f} nm".format(dl))

(x, y) = read_from_file("lum_blanche_1.txt")
(w, f, s) = fourierTransformInterferogram(x,y)
df = f[1] - f[0]
dl = 0.500*0.500*df*1000
#plotCombinedFigures(x, y, w, s, left=400, right=700, title="Lumière blanche 2, resolution {0:0.2f} nm".format(dl))

(x, y) = read_from_file("lum_blanche_2.txt")
(w, f, s) = fourierTransformInterferogram(x,y)
df = f[1] - f[0]
dl = 0.500*0.500*df*1000
#plotCombinedFigures(x, y, w, s, left=400, right=700, title="Lumière blanche 2, resolution {0:0.2f} nm".format(dl))

(x, y) = read_from_file("lum_blanche_3.txt")
(w, f, s) = fourierTransformInterferogram(x,y)
df = f[1] - f[0]
dl = 0.500*0.500*df*1000
#plotCombinedFigures(x, y, w, s, left=0, right=1000, title="Lumière blanche 3, resolution {0:0.2f} nm".format(dl))

(x, y) = read_from_file("lum_blanche_5.txt")
(w, f, s) = fourierTransformInterferogram(x,y)
df = f[1] - f[0]
dl = 0.500*0.500*df*1000
#plotCombinedFigures(x, y, w, s, left=0, right=10000, title="Lumière blanche 5, resolution {0:0.2f} nm".format(dl))

(x, y) = read_from_file("lum_blanche_6.txt")
(w, f, s) = fourierTransformInterferogram(x,y)
df = f[1] - f[0]
dl = 0.500*0.500*df*1000
plotCombinedFigures(x, y, w, s, left=0, right=1300, title="Lumière blanche, resolution {0:0.2f} nm".format(dl))

(x, y) = read_from_file("lum_sodium_1.txt")
(w, f, s) = fourierTransformInterferogram(x,y)
df = f[1] - f[0]
dl = 0.500*0.500*df*1000
plotCombinedFigures(x, y, w, s, left=3900, right=10000, title="Lumière sodium, resolution {0:0.2f} nm".format(dl))

(x, y) = read_from_file("lum_mercure_1.txt")
(w, f, s) = fourierTransformInterferogram(x,y)
df = f[1] - f[0]
dl = 0.500*0.500*df*1000
#plotCombinedFigures(x, y, w, s, left=100, right=800, title="Lumière Hg 1, resolution {0:0.2f} nm".format(dl))

(x, y) = read_from_file("lum_mercure_2.txt")
(w, f, s) = fourierTransformInterferogram(x,y)
df = f[1] - f[0]
dl = 0.500*0.500*df*1000
plotCombinedFigures(x, y, w, s, left=100, right=800, title="Lumière au mercure, resolution {0:0.2f} nm".format(dl))

(x, y) = read_from_file("lum_jaune_0.txt")
(w, f, s) = fourierTransformInterferogram(x,y)
df = f[1] - f[0]
dl = 0.500*0.500*df*1000
#plotCombinedFigures(x, y, w, s, left=400, right=7000, title="Lumière jaune, resolution {0:0.2f} nm".format(dl))

(x, y) = read_from_file("autre_jaune.txt")
(w, f, s) = fourierTransformInterferogram(x,y)
df = f[1] - f[0]
dl = 0.500*0.500*df*1000
plotCombinedFigures(x, y, w, s, left=100, right=800, title="Lumière jaune, resolution {0:0.2f} nm".format(dl))


plt.show() #J: À la fin pour afficher tous les graphiques en même temps