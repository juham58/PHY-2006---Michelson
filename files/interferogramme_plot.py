import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import *

""" Ce script fait la transformée de Fourier d'un fichier de données obtenu lors du laboratoire de l'Interféromètre
de Michelson fait dans le cours PHY-2006 - Optique expérimentale.

Le code est adapté de celui de Daniel Côté dispobile à l'adresse:
https://github.com/dccote/Enseignement/tree/master/ATELIER/FFT

J'ai changé la fonction d'extraction de données qui était utilisée auparavant pour la fonction genfromtxt de numpy qui
fonctionne mieux selon moi. J'ai ajouté un facteur 2 au dx de la fonction fourierTransformInterferogram pour prendre
en compte le trajet deux fois plus long de la lumière par rapport au déplacement du miroir M1.

J'ai aussi ajouté une petite fonction plot_ocean_optics qui permet de lire et d'afficher les données du spectrogramme
Ocean Optics utilisé en laboratoire comme référence.
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
	dx = 2*(x[1]-x[0]) # on obtient dx, on suppose equidistant J: ajout d'un facteur 2 pour le déplacement du miroir
	N = len(x)     # on obtient N directement des données
	frequencies = fftfreq(N, dx) # Cette fonction est fournie par numpy
	wavelengths = 1/frequencies-float(200e-9)  # Les fréquences en µm^-1 sont moins utiles que lambda en µm
	return (wavelengths, frequencies, spectrum)



def plotCombinedFigures(x, y, w, s, title="", left=400, right=800):
	""""
	On met l'interferogramme et le spectre sur la meme page.
	"""
	fig, (axes, axesFFT) = plt.subplots(2,1,figsize=(10, 10))
	axes.plot(x, y, '-m')
	axes.set_title("Interferogramme")
	axesFFT.plot(w*1000, abs(s), '.-m')
	axesFFT.set_xlim(left=left, right=right)
	axesFFT.set_xlabel("Longueur d'onde [nm]")
	axesFFT.set_title(title)



(x, y) = read_from_file("laser_he_ne_3.txt")
(w, f, s) = fourierTransformInterferogram(x,y)
df = f[1] - f[0]
dl = 0.500*0.500*df*1000
plotCombinedFigures(x, y, w, s, left=300, right=800, title="Laser He-Ne, resolution {0:0.2f} nm".format(dl))

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

(x, y) = read_from_file("lum_mercure_2.txt")
(w, f, s) = fourierTransformInterferogram(x,y)
df = f[1] - f[0]
dl = 0.500*0.500*df*1000
plotCombinedFigures(x, y, w, s, left=100, right=800, title="Lumière au mercure, resolution {0:0.2f} nm".format(dl))

(x, y) = read_from_file("autre_jaune.txt")
(w, f, s) = fourierTransformInterferogram(x,y)
df = f[1] - f[0]
dl = 0.500*0.500*df*1000
plotCombinedFigures(x, y, w, s, left=100, right=800, title="Lumière jaune, resolution {0:0.2f} nm".format(dl))


plt.show() #J: À la fin pour afficher tous les graphiques en même temps