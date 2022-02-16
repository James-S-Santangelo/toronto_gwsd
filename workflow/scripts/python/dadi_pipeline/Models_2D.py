import numpy
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum

'''
Models for testing one population scenari.
'''
def no_div (notused, ns, pts):
    """
    Standard neutral model, populations never diverge.
    """

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def no_div_bot(params, ns, pts):
    """
    Model with no population split, but the whole population goes
    through a bottleneck of depth nuB at time TB + TF followed by a
    recovery to size nuF at time TF
    """
    nuB, nuF, TB, TF = params
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)

    phi = Integration.one_pop(phi, xx, TB, nuB)
    phi = Integration.one_pop(phi, xx, TF, nuF)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def no_div_growth(params, ns, pts):
    """
    Model with no population split, but the whole population grows
    exponentionally to reach size nu at in the present
    """
    nu, T = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)

    nu_func = lambda t : numpy.exp(numpy.log(nu) * t/T)
    phi = Integration.one_pop(phi, xx, T, nu_func)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def no_div_bot_growth(params, ns, pts):
    """
    Model with no population split, but the whole population goes
    through a bottleneck of depth nuB at time TB + TF followed by an
    exponential growth to size nuF at time TF
    """
    nuB, nuG, TB, TG = params
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)

    phi = Integration.one_pop(phi, xx, TB, nuB)
    nu_func = lambda t : numpy.exp(numpy.log(nuB) * t/TG)
    phi = Integration.one_pop(phi, xx, TG, nu_func)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


'''
Models for testing two population scenarios.
'''

### Split no change in pop size ################################

def split_no_mig(params, ns, pts):
    """
    Split into two populations, no migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    """
    nu1, nu2, T = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def split_sym_mig(params, ns, pts):
    """
    Split into two populations, with symmetric migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    m: Migration rate between populations (2*Na*m)
    """
    nu1, nu2, m, T = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def split_asym_mig(params, ns, pts):
    """
    Split into two populations, with different migration rates.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
	"""
    nu1, nu2, m12, m21, T = params
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m12, m21=m21)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs    
################################################################

### Split bottleneck+ exp recovery in pop1 = URBAN ##################

def split_bot_urb_no_mig(params, ns, pts):
    """
    Split into two populations, no migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu1B: The bottleneck size for pop1
    nu1F: The final size for pop1
    Ts: The time between the split and bottleneck
    Tb: The scaled time between the bottleneck in pop 1 and present.
    """
    nu1, nu2, nu1B, nu1F, Tb, Ts= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    #split in 2 pop
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=0, m21=0)
    # function for bottleneck then exp recov in pop 1
    nu1_func = lambda t: nu1B*(nu1F/nu1B)**(t/Tb)
    phi = Integration.two_pops(phi, xx, Tb, nu1F=nu1_func, nu2=nu2, 
                                    m12=0, m21=0)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def split_bot_urb_sym_mig(params, ns, pts):
    """
    Split into two populations, with symetric migration

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu1B: The bottleneck size for pop1
    nu1F: The final size for pop1
    m: migration rate (constant)
    Ts: The time between the split and bottleneck
    Tb: The scaled time between the bottleneck in pop 1 and present.
    """
    nu1, nu2, nu1B, nu1F, m, Tb, Ts= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    #split in 2 pop
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=m, m21=m)
    # function for bottleneck then exp recov in pop 1
    nu1_func = lambda t: nu1B*(nu1F/nu1B)**(t/Tb)
    phi = Integration.two_pops(phi, xx, Tb, nu1F=nu1_func, nu2=nu2, 
                                    m12=m, m21=m)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def split_bot_urb_asym_mig(params, ns, pts):
    """
    Split into two populations, with symetric migration

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu1B: The bottleneck size for pop1
    nu1F: The final size for pop1
    m12: migration rate from pop 1 to 2 (constant)
    m21: migration rate from pop 2 to 1 (constant)
    Ts: The time between the split and bottleneck
    Tb: The scaled time between the bottleneck in pop 1 and present.
    """
    nu1, nu2, nu1B, nu1F, m12, m21, Tb, Ts= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    #split in 2 pop
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    # function for bottleneck then exp recov in pop 1
    nu1_func = lambda t: nu1B*(nu1F/nu1B)**(t/Tb)
    phi = Integration.two_pops(phi, xx, Tb, nu1F=nu1_func, nu2=nu2, 
                                    m12=m12, m21=m21)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs  
################################################################


### Split bottleneck+ exp recovery in pop2 = RURAL ##################

def split_bot_rur_no_mig(params, ns, pts):
    """
    Split into two populations, no migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu2B: The bottleneck size for pop2
    nu2F: The final size for pop2
    Ts: The time between the split and bottleneck
    Tb: The scaled time between the bottleneck in pop 2 and present.
    """
    nu1, nu2, nu2B, nu2F, Tb, Ts= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    #split in 2 pop
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=0, m21=0)
    # function for bottleneck then exp recov in pop 1
    nu2_func = lambda t: nu2B*(nu2F/nu2B)**(t/Tb)
    phi = Integration.two_pops(phi, xx, Tb, nu1=nu1, nu2F=nu2_func, 
                                    m12=0, m21=0)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def split_bot_rur_sym_mig(params, ns, pts):
    """
    Split into two populations, with symetric migration

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu2B: The bottleneck size for pop2
    nu2F: The final size for pop2
    m: migration rate (constant)
    Ts: The time between the split and bottleneck
    Tb: The scaled time between the bottleneck in pop 2 and present.
    """
    nu1, nu2, nu2B, nu2F, m, Tb, Ts= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    #split in 2 pop
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=m, m21=m)
    # function for bottleneck then exp recov in pop 1
    nu2_func = lambda t: nu2B*(nu2F/nu2B)**(t/Tb)
    phi = Integration.two_pops(phi, xx, Tb, nu1=nu1, nu2F=nu2_func, 
                                    m12=m, m21=m)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def split_bot_rur_asym_mig(params, ns, pts):
    """
    Split into two populations, with symetric migration

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu2B: The bottleneck size for pop1
    nu2F: The final size for pop1
    m12: migration rate from pop 1 to 2 (constant)
    m21: migration rate from pop 2 to 1 (constant)
    Ts: The time between the split and bottleneck
    Tb: The scaled time between the bottleneck in pop 2 and present.
    """
    nu1, nu2, nu2B, nu2F, m12, m21, Tb, Ts= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    #split in 2 pop
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    # function for bottleneck then exp recov in pop 1
    nu2_func = lambda t: nu2B*(nu2F/nu2B)**(t/Tb)
    phi = Integration.two_pops(phi, xx, Tb, nu1=nu1, nu2F=nu2_func, 
                                    m12=m12, m21=m21)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs  
################################################################
'''

# SIMILAR TO THE PREVIOUS MODELS ###
### Split and exp growth in pop1 = URBAN #######################

def split_growth_urb_no_mig(params, ns, pts):
    """
    Split into two populations, no migration. 
    Pop1 experiencing exponential growth right after split.
    s: proportion of Nref contributing to population 1 after split.
    nu1: Final size of population 1 after growth(present size).
    nu2: Size of population 2 after split and now (constant).
    T: The time between the split (=start of exp growth for pop1) and present
    """
    nu1, nu2, s, T= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    #split in 2 pop
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    # function for exp growth in pop 1 starting just at the split T
    nu1_func = lambda t: s*(nu1/s)**(t/T)
    phi = Integration.two_pops(phi, xx, T, nu1=nu1_func, nu2=nu2)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def split_growth_urb_sym_mig(params, ns, pts):
    """
    Split into two populations, symetric migration. 
    Pop1 experiencing exponential growth right after split.
    nu1: Size of population 1 after split.
    nu1F: Final size of population 1 after growth(present size).
    nu2: Size of population 2 after split and now (constant).
    m: migration rate (constant)
    T: The time between the split (=start of exp growth for pop1) and present
    """
    nu1, nu2, nu1F, m, T= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    #split in 2 pop
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=m, m21=m)
    # function for exp growth in pop 1 starting just at the split T
    nu1_func = lambda t: nu1*(nu1F/nu1)**(t/T)
    phi = Integration.two_pops(phi, xx, T, nu1F=nu1_func, nu2=nu2, m12=m, m21=m)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def split_growth_urb_asym_mig(params, ns, pts):
    """
    Split into two populations, asymetric migration. 
    Pop1 experiencing exponential growth right after split.
    nu1: Size of population 1 after split.
    nu1F: Final size of population 1 after growth(present size).
    nu2: Size of population 2 after split and now (constant).
    m12: migration rate from pop1 to pop2 (constant).
    m21: migration rate from pop2 to pop1 (constant).
    T: The time between the split (=start of exp growth for pop1) and present
    """
    nu1, nu2, nu1F, m12, m21, T= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    #split in 2 pop
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m12, m21=m21)
    # function for exp growth in pop 1 starting just at the split T
    nu1_func = lambda t: nu1*(nu1F/nu1)**(t/T)
    phi = Integration.two_pops(phi, xx, T, nu1F=nu1_func, nu2=nu2, m12=m12, m21=m21)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

################################################################

### Split and exp growth in pop2 = RURAL #######################

def split_growth_rur_no_mig(params, ns, pts):
    """
    Split into two populations, no migration. 
    Pop2 experiencing exponential growth right after split.
    nu1: Size of population 1 after split and now (constant).
    nu2: Size of population 2 after split.
    nu2F: Final size of population 2 (present size)
    
    T: The time between the split (=start of exp growth for pop1) and present
    """
    nu1, nu2, nu2F, T= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    #split in 2 pop
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=0, m21=0)
    # function for exp growth in pop 1 starting just at the split T
    nu2_func = lambda t: nu2*(nu2F/nu2)**(t/T)
    phi = Integration.two_pops(phi, xx, T, nu1, nu2F=nu2_func, m12=0, m21=0)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def split_growth_rur_sym_mig(params, ns, pts):
    """
    Split into two populations, symetric migration. 
    Pop2 experiencing exponential growth right after split.
    nu1: Size of population 1 after split and now (constant).
    nu2: Size of population 2 after split.
    nu2F: Final size of population 2 (present size).
    m: migration rate (constant)
    T: The time between the split (=start of exp growth for pop1) and present
    """
    nu1, nu2, nu2F, m, T= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    #split in 2 pop
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m, m21=m)
    # function for exp growth in pop 1 starting just at the split T
    nu2_func = lambda t: nu2*(nu2F/nu2)**(t/T)
    phi = Integration.two_pops(phi, xx, T, nu1, nu2F=nu2_func, m12=m, m21=m)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def split_growth_rur_asym_mig(params, ns, pts):
    """
    Split into two populations, asymetric migration. 
    Pop2 experiencing exponential growth right after split.
    nu1: Size of population 1 after split and now (constant).
    nu2: Size of population 2 after split.
    nu2F: Final size of population 2 (present size).
    m12: migration rate from pop1 to pop2 (constant)
    m21: migration rate from pop2 to pop1
    T: The time between the split (=start of exp growth for pop1) and present
    """
    nu1, nu2, nuF, m12, m21, T= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    #split in 2 pop
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m12, m21=m21)
    # function for exp growth in pop 1 starting just at the split T
    nu2_func = lambda t: nu2*(nu2F/nu2)**(t/T)
    phi = Integration.two_pops(phi, xx, T, nu1, nu2F=nu2_func, m12=m12, m21=m21)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

################################################################

'''





''' WORK IN PROGRESS
### Split bottleneck+ exp recovery in pop1 and pop2 #############

def split_bot_urb_rur_no_mig(params, ns, pts):
    """
    Split into two populations, no migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu1B: The bottleneck size for pop1
    nu2B: The bottleneck size for pop1
    nu2B: The bottleneck size for pop2
    nu2F: The final size for pop2
    Ts: The time between the split and bottleneck for pop 1
    T1b: The scaled time between the bottleneck in pop 1 and present.
    T2b: The scaled time between the bottleneck in pop 2 and present.
    """
    nu1, nu2, nu1B, nu1F, nu2B, nu2F, T1b, T2b, Ts= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    #split in 2 pop
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=0, m21=0)
    # function for bottleneck then exp recov in pop 1
    nu1_func = lambda t: nu1B*(nu1F/nu1B)**(t/T1b)
    nu2_func = lambda t: nu2B*(nu2F/nu2B)**(t/T2b)
    phi = Integration.two_pops(phi, xx, T1b, nu1=nu1_func, nu2=nu2, 
                                    m12=0, m21=0)
    phi = Integration.two_pops(phi, xx, T2b, nu1=nu1, nu2=nu2_func, 
                                    m12=0, m21=0)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def split_bot_rur_sym_mig(params, ns, pts):
    """
    Split into two populations, no migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu1B: The bottleneck size for pop1
    nu2B: The bottleneck size for pop1
    nu2B: The bottleneck size for pop2
    nu2F: The final size for pop2
    m: migration rate
    Ts: The time between the split and bottleneck for pop 1
    T1b: The scaled time between the bottleneck in pop 1 and present.
    T2b: The scaled time between the bottleneck in pop 2 and present.
    """
    nu1, nu2, nu1B, nu1F, nu2B, nu2F, T1b, T2b, Ts= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    #split in 2 pop
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=m, m21=m)
    # function for bottleneck then exp recov in pop 1
    nu1_func = lambda t: nu1B*(nu1F/nu1B)**(t/T1b)
    nu2_func = lambda t: nu2B*(nu2F/nu2B)**(t/T2b)
    phi = Integration.two_pops(phi, xx, T1b, nu1=nu1_func, nu2=nu2, 
                                    m12=m, m21=m)
    phi = Integration.two_pops(phi, xx, T2b, nu1=nu1, nu2=nu2_func, 
                                    m12=m, m21=m)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def split_bot_rur_asym_mig(params, ns, pts):
    """
    Split into two populations, no migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu1B: The bottleneck size for pop1
    nu2B: The bottleneck size for pop1
    nu2B: The bottleneck size for pop2
    nu2F: The final size for pop2
    m12: migration rate from pop 1 to 2 (constant)
    m21: migration rate from pop 2 to 1 (constant)
    Ts: The time between the split and bottleneck for pop 1
    T1b: The scaled time between the bottleneck in pop 1 and present.
    T2b: The scaled time between the bottleneck in pop 2 and present.
    """
    nu1, nu2, nu1B, nu1F, nu2B, nu2F, T1b, T2b, Ts= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    #split in 2 pop
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    # function for bottleneck then exp recov in pop 1
    nu1_func = lambda t: nu1B*(nu1F/nu1B)**(t/T1b)
    nu2_func = lambda t: nu2B*(nu2F/nu2B)**(t/T2b)
    phi = Integration.two_pops(phi, xx, T1b, nu1=nu1_func, nu2=nu2, 
                                    m12=m12, m21=m21)
    phi = Integration.two_pops(phi, xx, T2b, nu1=nu1, nu2=nu2_func, 
                                    m12=m12, m21=m21)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs 
################################################################
'''
