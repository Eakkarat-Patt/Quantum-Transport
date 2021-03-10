import numpy as np
import matplotlib.pyplot as plt
plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
plt.rc('text',usetex=True)


phi = np.linspace(-np.pi / 2, np.pi / 2, 10000)
Volt = np.linspace(170e-3,300e-3,10000)
hb = 6.58212e-16
vF = 1e6
L = 1e-7

"""
Only region in the middle is tilted
"""
def sign1(Ef):
    s = np.sign(Ef)
    return s


def sign2(Ef, U):
    sp = np.sign(Ef - U)
    return sp


def sign3(B0):
    eta = np.sign(B0)
    return eta


def wave_vector_ky(Ef, B0, phi):
    ky = (Ef / (hb * vF) * np.sin(phi)) + (sign3(B0) * np.sqrt((abs(B0)) / hb))
    return ky

def wave_vector_k(Ef, e, phi):
    k = Ef / sign1(Ef) * hb * vF * (1-e * np.sin(phi))
    return k


def wave_vector_q(Ef, U, B0, e, phi):
    q = ((Ef - U) / (hb * vF * sign2(Ef, U))) + ((wave_vector_ky(Ef, B0, phi) * e) / (sign2(Ef, U))+ 0j)
    return q


def wave_vector_qx(Ef, U, B0, e, phi):
    qx = np.sqrt(wave_vector_q(Ef, U, B0, e, phi) ** 2 - (wave_vector_ky(Ef, B0, phi)) ** 2 + 0j)
    return qx


def refrac(Ef, U, B0, e, phi):
    theta = np.arcsin(wave_vector_ky(Ef, B0, phi) / wave_vector_q(Ef, U, B0, e, phi)+ 0j)
    return theta


def trans(Ef, U, B0, e, phi):  # transmission probability
    tp1 = (np.cos(refrac(Ef, U, B0, e, phi)) * np.cos(phi)) ** 2
    tp2 = (np.cos(L * wave_vector_qx(Ef, U, B0, e, phi)) * np.cos(refrac(Ef, U, B0, e, phi)) * np.cos(phi)) ** 2
    tp3 = np.sin(L * wave_vector_qx(Ef, U, B0, e, phi)) ** 2 * (1 - (sign1(Ef) * sign2(Ef, U) * np.sin(refrac(Ef, U, B0, e, phi)) * np.sin(phi))) ** 2
    TP = tp1 / (tp2 + tp3)
    return TP


"""
All three regions are tilted
"""

def wave_vector_k2(Ef, e, phi):
    k2 = Ef/(hb * vF * (1 + e * np.cos((np.pi/2) + phi)))
    return k2

def wave_vector_q2(Ef, U, e, phi):
    q2 = ((Ef-U)/(hb * vF * sign2(Ef, U))) + ((wave_vector_k2(Ef, e, phi) * np.sin(phi) * e)/(sign2(Ef, U)))
    return q2

def refrac2(Ef, U , e, phi):
    theta2 = np.arcsin(wave_vector_k2(Ef, e, phi) * np.sin(phi)/(wave_vector_q2(Ef, U, e, phi)) + 0j)
    return theta2

def wave_vector_qx2(Ef, U, e, phi):
    qx2 = np.sqrt(wave_vector_q2(Ef, U, e, phi)**2 -(wave_vector_k2(Ef, e, phi)*np.sin(phi))**2 + 0j)
    return qx2

def trans2(Ef, U, e, phi):
    tp4 = (np.cos(refrac2(Ef, U, e, phi)) * np.cos(phi))**2
    tp5 = (np.cos(L * wave_vector_qx2(Ef, U, e, phi)) * np.cos(refrac2(Ef, U, e, phi)) * np.cos(phi))**2
    tp6 = np.sin(L * wave_vector_qx2(Ef, U, e, phi))**2 * (1 - (sign1(Ef) * sign2(Ef, U) * np.sin(refrac2(Ef, U, e, phi)) * np.sin(phi)))**2
    TP2 =  tp4/(tp5+tp6)
    return TP2



def tp_plot(Ef, U, e1, e2, e3, l, h, pr):
    """
    Polar plot: transmission prob vs angle
    Study the transmission shift
    set l = h = 8
    """
    fig = plt.figure(figsize= (l, h))
    ax = plt.subplot(111, projection="polar")
    ax.set_rlabel_position(90)
    ax.plot(phi, trans(Ef, U, 0, e1, phi), linewidth="3", linestyle="-")
    ax.plot(phi, trans(Ef, U, 0, e2, phi), linewidth="3", linestyle="--")
    ax.plot(phi, trans(Ef, U, 0, e3, phi), linewidth ="3",linestyle = "-.")
    #ax.legend(("$w_0 = {}$ ".format(e1),"$w_0 = {}$".format(e2), "$w_0 = {}$".format(e3)), bbox_to_anchor=(0.7, 1), frameon = False, labelspacing = 1, fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_title("transmission probability Ef = {} U = {} e1 = {} e2 = {} e3 = {}".format(Ef, U, e1, e2, e3))
    ax.set_rmax(1)
    ax.set_thetamin(-90)
    ax.set_thetamax(90)
    if pr ==1:
        plt.savefig("/Users/eakkarat/Google Drive/Thesis (Master)/Thesis/fig/tunneling shift/U{}.png".format(U), dpi=1000)
    elif pr ==0:
        plt.show()

def tp_plot3(Ef, U, B, w, l, h,lx,ly, pr):
    """
    Polar plot: transmission prob vs angle
    Compare the tp: real vs pseudo B field
    set l = h = 8
    """
    fig = plt.figure(figsize= (l, h))
    ax = plt.subplot(111, projection="polar")
    ax.set_rlabel_position(90)
    ax.plot(phi, trans(Ef, U, 0, w, phi), linewidth="3", linestyle="-")
    ax.plot(phi, trans(Ef, U, B, 0, phi), linewidth="3", linestyle="--")
    ax.legend(("$w_0 = {}$ ".format(w),"$B = {}$".format(B)), bbox_to_anchor=(lx, ly), frameon = False, labelspacing = 1, fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_title("transmission probability Ef = {} U = {} B = {} w = {}".format(Ef, U, B, w))
    ax.set_rmax(1)
    ax.set_thetamin(-90)
    ax.set_thetamax(90)
    if pr ==1:
        plt.savefig("/Users/eakkarat/Google Drive/Thesis (Master)/Thesis/fig/pseudo B field/Ef{} U{} B{} w{}.png".format(Ef, U, B, w), dpi=1000)
    elif pr ==0:
        plt.show()

"""Cartesian plot: transmission prob vs angle"""
def tp_plot_car(Ef, U, e1, e2, e3, l, h, pr):  # cartesian
    """Cartesian plot: transmission prob vs angle"""

    fig = plt.figure(figsize=(l, h))
    ax = plt.subplot(111)
    ax.plot(phi, trans(Ef, U, 0, e1, phi), linewidth="3", linestyle="-")
    ax.plot(phi, trans(Ef, U, 0, e2, phi), linewidth="3", linestyle="--")
    ax.plot(phi, trans(Ef, U, 0, e3, phi), linewidth="3", linestyle="-.")
    #ax.legend((" ", " ", " "), bbox_to_anchor=(lx, ly), frameon=False, labelspacing = 1.5, fontsize='x-large')
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_title(
        "transmission probability Ef = {} U = {} e1 = {} e2 = {} e3 = {}".format(Ef, U, e1, e2, e3))
    #ax.set_xlim([0.17, 0.3])
    #ax.set_ylim([0, 1.1])
    if pr ==1:
        plt.savefig("./fig/tp Ef = {} B = {} e1 = {} e2 = {} e3 = {}.png".format(Ef, B, e1, e2, e3), dpi=1000)
    elif pr == 0:
        plt.show()

def tp_plot_car2(Ef, e1, e2, e3, l, h, lx, ly, pr):
    """Cartesian plot at fixed angle
    l = 13, h = 9
    """
    fig = plt.figure(figsize=(l,h))
    ax = plt.subplot(111)
    ax.plot(Volt, trans(Ef, Volt, 0, e1, np.pi/4), linewidth="3", linestyle="-")
    ax.plot(Volt, trans(Ef, Volt, 0, e2, np.pi/4), linewidth="3", linestyle="--")
    ax.plot(Volt, trans(Ef, Volt, 0, e3, np.pi/4), linewidth="3", linestyle="-.")
    ax.legend(("$w_0 = {}$".format(e1), "$w_0 = {}$".format(e2), "$w_0 = {}$".format(e3)), bbox_to_anchor=(lx, ly), frameon=False, labelspacing = 1.0, fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_title(
        "transmission probability Ef = {} U = {} e1 = {} e2 = {} e3 = {}".format(Ef, Volt, e1, e2, e3))
    fig.set_size_inches(l, h)
    ax.set_xlabel("$U$ (eV)", fontsize = 20)
    ax.set_ylabel("$T$", fontsize = 20)
    ax.set_xlim([0.17, 0.3])
    ax.set_ylim([0, 1.1])
    if pr ==1:
        plt.savefig("./fig/tp Ef = {} e1 = {} e2 = {} e3 = {}.png".format(Ef, e1, e2, e3), dpi=1000)
    elif pr == 0:
        plt.show()
"""
Contour plot of transmission (phi vs ec)
"""
def tp_plot_contour(Ef, U, l, h, pr):
    """l = 8, h =7
    """

    ec = np.linspace(0, 1, 1000)
    Phi, Ec = np.meshgrid(phi, ec)
    levels = np.arange(0, 1, 0.005)
    fig = plt.figure(figsize=(l, h))
    ax = plt.subplot(111)
    ax.contourf(Phi, Ec, trans(Ef, U,0, Ec, Phi), levels, colors = 'k')
    ax.set_title("tp contour Ef = {} U = {}".format(Ef, U))
    contour_filled = plt.contourf(Phi, Ec, trans(Ef, U,0, Ec, Phi), levels)
    ax.tick_params(labelbottom=True, labeltop=True, labelleft=True, labelright=True,
                   bottom=True, top=True, left=True, right=True, labelsize = 15)
    cbar = plt.colorbar(contour_filled, pad = 0.1)
    cbar.ax.tick_params(labelsize = 15)
    if pr == 1:
        plt.savefig("./fig/tp/contour Ef = {} U = {}.png".format(Ef, U),
                dpi=1000)
    elif pr == 0:
        plt.show()


def tp_plot2(Ef, U, e1, l, h):
    """"
     Polar plot: transmission prob vs angle with and without mismatch effect
    set l = h = 8
    """
    fig = plt.figure(figsize=(l, h))
    ax = plt.subplot(111, projection="polar")
    ax.set_rlabel_position(90)
    ax.plot(phi, trans(Ef, U, 0, e1, phi), linewidth="3", linestyle="-")
    ax.plot(phi, trans2(Ef, U, e1, phi), linewidth="3", linestyle="--")
    #ax.legend((" ", " ", " "), bbox_to_anchor=(0.35, 0.5), frameon=False, labelspacing=1, fontsize='x-large')
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_title("transmission probability Ef = {} U = {} e1 = {}".format(Ef, U, e1))
    ax.set_rmax(1)
    ax.set_thetamin(-90)
    ax.set_thetamax(90)
    plt.savefig("./fig/tp/tp Ef = {} U = {} e1 = {}.png".format(Ef, U, e1), dpi=1000)
    plt.show()