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
    phi = Integration.two_pops(phi, xx, Tb, nu1=nu1_func, nu2=nu2, 
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
    phi = Integration.two_pops(phi, xx, Tb, nu1=nu1_func, nu2=nu2, 
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
    phi = Integration.two_pops(phi, xx, Tb, nu1=nu1_func, nu2=nu2, 
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
    phi = Integration.two_pops(phi, xx, Tb, nu1=nu1, nu2=nu2_func, 
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
    phi = Integration.two_pops(phi, xx, Tb, nu1=nu1, nu2=nu2_func, 
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
    phi = Integration.two_pops(phi, xx, Tb, nu1=nu1, nu2=nu2_func, 
                                    m12=m12, m21=m21)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs  

################################################################


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
    nu1, s, T= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    # function for exp growth in pop 1 starting just at the split T
    nu1_func = lambda t: s*(nu1/s)**(t/T)
    # pop1 immidiately starts to exponentialy grow after split
    phi = Integration.two_pops(phi, xx, T, nu1=nu1_func, nu2=1-s)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def split_growth_urb_sym_mig(params, ns, pts):
    """
    Split into two populations, symmetric migration. 
    Pop1 experiencing exponential growth right after split.
    s: proportion of Nref contributing to population 1 after split.
    nu1: Final size of population 1 after growth(present size).
    nu2: Size of population 2 after split and now (constant).
    T: The time between the split (=start of exp growth for pop1) and presentt).
    m: migration rate (constant)
    """

    nu1, s, T, m = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    # function for exp growth in pop 1 starting just at the split T
    nu1_func = lambda t: s*(nu1/s)**(t/T)
    # pop1 immidiately starts to exponentialy grow after split
    phi = Integration.two_pops(phi, xx, T, nu1=nu1_func, nu2=1-s, m12=m, m21=m)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def split_growth_urb_asym_mig(params, ns, pts):
    """
    Split into two populations, symmetric migration. 
    Pop1 experiencing exponential growth right after split.
    s: proportion of Nref contributing to population 1 after split.
    nu1: Final size of population 1 after growth(present size).
    nu2: Size of population 2 after split and now (constant).
    m12: migration rate from pop1 to pop2 (constant).
    m21: migration rate from pop2 to pop1 (constant).
    T: The time between the split (=start of exp growth for pop1) and present
    """
    nu1, s, T, m12, m21 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    # function for exp growth in pop 1 starting just at the split T
    nu1_func = lambda t: s*(nu1/s)**(t/T)
    # pop1 immidiately starts to exponentialy grow after split
    phi = Integration.two_pops(phi, xx, T, nu1=nu1_func, nu2=1-s, m12=m12, m21=m21)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs
    
################################################################

### Split and exp growth in pop2 = RURAL #######################

def split_growth_rur_no_mig(params, ns, pts):
    """
    Split into two populations, no migration. 
    Pop2 experiencing exponential growth right after split.
    s: proportion of Nref contributing to population 2 after split.
    nu1: Size of population 1 after split (constant = 1-s).
    nu2: Size of population 2 after split and growth (present size).
    T: The time between the split (=start of exp growth for pop2) and present
    """
    nu2, s, T= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    # function for exp growth in pop 2 starting just at the split T
    nu2_func = lambda t: s*(nu2/s)**(t/T)
    # pop1 immidiately starts to exponentialy grow after split
    phi = Integration.two_pops(phi, xx, T, nu1=1-s, nu2=nu2_func)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs



def split_growth_rur_sym_mig(params, ns, pts):
    """
    Split into two populations, symmetric migration. 
    Pop2 experiencing exponential growth right after split.
    s: proportion of Nref contributing to population 1 after split.
    nu1: Size of population 1 after split (constant = 1-s).
    nu2: Size of population 2 after split and growth (present size).
    T: The time between the split (=start of exp growth for pop2) and presentt).
    m: migration rate (constant)
    """
    nu2, s, T, m = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    # function for exp growth in pop 2 starting just at the split T
    nu2_func = lambda t: s*(nu2/s)**(t/T)
    # pop1 immidiately starts to exponentialy grow after split
    phi = Integration.two_pops(phi, xx, T, nu1=1-s, nu2=nu2_func, m12=m, m21=m)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def split_growth_rur_asym_mig(params, ns, pts):
    """
    Split into two populations, asymmetric migration. 
    Pop2 experiencing exponential growth right after split.
    nu1: Size of population 1 after split (constant = 1-s).
    nu2: Size of population 2 after split and growth (present size).
    m12: migration rate from pop1 to pop2 (constant).
    m21: migration rate from pop2 to pop1 (constant).
    T: The time between the split (=start of exp growth for pop2) and present
    """
    nu2, s, T, m12, m21= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    # function for exp growth in pop 2 starting just at the split T
    nu2_func = lambda t: s*(nu2/s)**(t/T)
    # pop1 immidiately starts to exponentialy grow after split
    phi = Integration.two_pops(phi, xx, T, nu1=1-s, nu2=nu2_func, m12=m12, m21=m21)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs
  
################################################################


### Split and immediate exp growth in pop1 = URBAN and then growth of pop2 #######################

def split_growth_urb_then_late_growth_rur_no_mig(params, ns, pts):
    """
    Split into two populations, no migration. 
    Pop1 experiencing exponential growth right after split.
    Pop2 starts its exponential growth later.
    s: proportion of Nref contributing to population 1 after split.
    nu1: Final size of population 1 after growth(present size).
    nu2: Size of population 2 after growth.
    Ts: The time between the split (=start of exp growth for pop1) and exp groth of pop2
    Tg2: time since the beggining of pop2 exp growth
    """
    nu1, nu2, s, Ts, Tg2= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    # function for exp growth in pop 1 starting just at the split T
    nu1_func = lambda t: s*(nu1/s)**(t/(Ts+Tg2))
    # pop1 immidiately starts to exponentialy grow after split (Ts+Tg2 unit of time ago)
    phi = Integration.two_pops(phi, xx, Ts, nu1=nu1_func, nu2=1-s)
    # Then, Tg2 units of time ago, pop2 starts growing exponentialy 
    nu2_func = lambda t: (1-s)*(nu2/(1-s))**(t/Tg2)
    phi = Integration.two_pops(phi, xx, Tg2, nu1=nu1_func, nu2=nu2_func)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def split_growth_urb_then_late_growth_rur_sym_mig(params, ns, pts):
    """
    Split into two populations, symmetric migration. 
    Pop1 experiencing exponential growth right after split.
    Pop2 starts its exponential growth later.
    s: proportion of Nref contributing to population 1 after split.
    nu1: Final size of population 1 after growth(present size).
    nu2: Size of population 2 after growth.
    Ts: The time between the split (=start of exp growth for pop1) and exp groth of pop2
    Tg2: time since the beggining of pop2 exp growth
    m: migration rate (constant)
    """

    nu1, nu2, s, Ts, Tg2, m= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    # function for exp growth in pop 1 starting just at the split T
    nu1_func = lambda t: s*(nu1/s)**(t/(Ts+Tg2))
    # pop1 immidiately starts to exponentialy grow after split (Ts+Tg2 unit of time ago)
    phi = Integration.two_pops(phi, xx, Ts, nu1=nu1_func, nu2=1-s, m12=m, m21=m)
    # Then, Tg2 units of time ago, pop2 starts growing exponentialy 
    nu2_func = lambda t: (1-s)*(nu2/(1-s))**(t/Tg2)
    phi = Integration.two_pops(phi, xx, Tg2, nu1=nu1_func, nu2=nu2_func, m12=m, m21=m)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def split_growth_urb_then_late_growth_rur_asym_mig(params, ns, pts):
    """
    Split into two populations, symmetric migration. 
    Pop1 experiencing exponential growth right after split.
    Pop2 starts its exponential growth later.
    s: proportion of Nref contributing to population 1 after split.
    nu1: Final size of population 1 after growth(present size).
    nu2: Size of population 2 after growth.
    Ts: The time between the split (=start of exp growth for pop1) and exp groth of pop2
    Tg2: time since the beggining of pop2 exp growth
    m12: migration rate from pop1 to pop2 (constant).
    m21: migration rate from pop2 to pop1 (constant).
    """

    nu1, nu2, s, Ts, Tg2, m12, m21= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    # function for exp growth in pop 1 starting just at the split T
    nu1_func = lambda t: s*(nu1/s)**(t/(Ts+Tg2))
    # pop1 immidiately starts to exponentialy grow after split (Ts+Tg2 unit of time ago)
    phi = Integration.two_pops(phi, xx, Ts, nu1=nu1_func, nu2=1-s, m12=m12, m21=m21)
    # Then, Tg2 units of time ago, pop2 starts growing exponentialy 
    nu2_func = lambda t: (1-s)*(nu2/(1-s))**(t/Tg2)
    phi = Integration.two_pops(phi, xx, Tg2, nu1=nu1_func, nu2=nu2_func, m12=m12, m21=m21)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs
################################################################

### Split and immediate exp growth in pop2 = RURAL and then growth of pop1 #######################

def split_growth_rur_then_late_growth_urb_no_mig(params, ns, pts):
    """
    Split into two populations, no migration. 
    Pop2 experiencing exponential growth right after split (from size s).
    Pop1 starts its exponential growth later.
    s: proportion of Nref contributing to population 2 after split.
    nu1: Final size of population 1 after growth(present size).
    nu2: Size of population 2 after growth.
    Ts: The time between the split (=start of exp growth for pop1) and exp groth of pop2
    Tg1: time since the beggining of pop1 exp growth
    """
    nu1, nu2, s, Ts, Tg1= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    # function for exp growth in pop 2 starting just at the split T
    nu2_func = lambda t: s*(nu2/s)**(t/(Ts+Tg1))
    # pop1 immidiately starts to exponentialy grow after split (Ts+Tg1 unit of time ago)
    phi = Integration.two_pops(phi, xx, Ts, nu1=1-s, nu2=nu2_func)
    # Then, Tg1 units of time ago, pop1 starts growing exponentialy 
    nu1_func = lambda t: (1-s)*(nu1/(1-s))**(t/Tg1)
    phi = Integration.two_pops(phi, xx, Tg1, nu1=nu1_func, nu2=nu2_func)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def split_growth_rur_then_late_growth_urb_sym_mig(params, ns, pts):
    """
    Split into two populations, symmetrci migration. 
    Pop2 experiencing exponential growth right after split (from size s).
    Pop1 starts its exponential growth later.
    s: proportion of Nref contributing to population 2 after split.
    nu1: Final size of population 1 after growth(present size).
    nu2: Size of population 2 after growth.
    Ts: The time between the split (=start of exp growth for pop1) and exp groth of pop2
    Tg1: time since the beggining of pop1 exp growth
    m: migration rate (constant)
    """
    nu1, nu2, s, Ts, Tg1, m= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
   # function for exp growth in pop 2 starting just at the split T
    nu2_func = lambda t: s*(nu2/s)**(t/(Ts+Tg1))
    # pop1 immidiately starts to exponentialy grow after split (Ts+Tg1 unit of time ago)
    phi = Integration.two_pops(phi, xx, Ts, nu1=1-s, nu2=nu2_func, m12=m, m21=m)
    # Then, Tg1 units of time ago, pop1 starts growing exponentialy 
    nu1_func = lambda t: (1-s)*(nu1/(1-s))**(t/Tg1)
    phi = Integration.two_pops(phi, xx, Tg1, nu1=nu1_func, nu2=nu2_func,m12=m, m21=m)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def split_growth_rur_then_late_growth_urb_asym_mig(params, ns, pts):
    """
    Split into two populations, symmetrci migration. 
    Pop2 experiencing exponential growth right after split (from size s).
    Pop1 starts its exponential growth later.
    s: proportion of Nref contributing to population 2 after split.
    nu1: Final size of population 1 after growth(present size).
    nu2: Size of population 2 after growth.
    Ts: The time between the split (=start of exp growth for pop1) and exp groth of pop2
    Tg1: time since the beggining of pop1 exp growth.
    m12: migration rate from pop1 to pop2 (constant).
    m21: migration rate from pop2 to pop1 (constant).
    """
    nu1, nu2, s, Ts, Tg1, m12, m21= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
   # function for exp growth in pop 2 starting just at the split T
    nu2_func = lambda t: s*(nu2/s)**(t/(Ts+Tg1))
    # pop1 immidiately starts to exponentialy grow after split (Ts+Tg1 unit of time ago)
    phi = Integration.two_pops(phi, xx, Ts, nu1=1-s, nu2=nu2_func, m12=m12, m21=m21)
    # Then, Tg1 units of time ago, pop1 starts growing exponentialy 
    nu1_func = lambda t: (1-s)*(nu1/(1-s))**(t/Tg1)
    phi = Integration.two_pops(phi, xx, Tg1, nu1=nu1_func, nu2=nu2_func, m12=m12, m21=m21)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs
################################################################

### Split and immediate exp growth in the two pops #############
def split_growth_urb_rur_no_mig(params, ns, pts):
    """
    Split into two populations, symmetrci migration. 
    Pop1 and Pop2  experiencing exponential growth right after split (from size s and 1-s).
    s: proportion of Nref contributing to population 1 after split.
    nu1: Size of population 1 after growth.
    nu2: Size of population 2 after growth.
    T: The time between the split and present
    """
    s, nu1, nu2, T = params
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t : s * (nu1/s) ** (t/T)
    nu2_func = lambda t : (1-s) * (nu2/(1-s)) ** (t/T)

    phi = Integration.two_pops(phi, xx, T, nu1_func, nu2_func)

    fs = Spectrum.from_phi(phi, ns, (xx, xx))
    return fs

def split_growth_urb_rur_sym_mig(params, ns, pts):
    """
    Split into two populations, symmetrci migration. 
    Pop1 and Pop2  experiencing exponential growth right after split (from size s and 1-s).
    s: proportion of Nref contributing to population 1 after split.
    nu1: Size of population 1 after growth.
    nu2: Size of population 2 after growth.
    T: The time between the split and present
    m: migration rate.
    """
    s, nu1, nu2, T, m = params
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t : s * (nu1/s) ** (t/T)
    nu2_func = lambda t : (1-s) * (nu2/(1-s)) ** (t/T)

    phi = Integration.two_pops(phi, xx, T, nu1_func, nu2_func, m12 = m, m21 = m)

    fs = Spectrum.from_phi(phi, ns, (xx, xx))
    return fs

def split_growth_urb_rur_asym_mig(params, ns, pts):
    """
    Split into two populations, symmetrci migration. 
    Pop1 and Pop2  experiencing exponential growth right after split (from size s and 1-s).
    s: proportion of Nref contributing to population 1 after split.
    nu1: Size of population 1 after growth.
    nu2: Size of population 2 after growth.
    T: The time between the split and present
    m12: migration rate from pop1 to pop2 (constant).
    m21: migration rate from pop2 to pop1 (constant).
    """
    s, nu1, nu2, T, m12, m21 = params
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t : s * (nu1/s) ** (t/T)
    nu2_func = lambda t : (1-s) * (nu2/(1-s)) ** (t/T)

    phi = Integration.two_pops(phi, xx, T, nu1_func, nu2_func, m12 = m12, m21 = m21)

    fs = Spectrum.from_phi(phi, ns, (xx, xx))
    return fs

################################################################

### Split, cst size then exp growth in pop2 the in pop1 ########

def split_bot_growth_rur_then_late_growth_urb_no_mig(params, ns, pts):
    """
    Split into two populations, no migration. 
    Pop2 experiencing exponential growth few moments after split (from size 1-s).
    Pop2 starts its exponential growth later.
    s: proportion of Nref contributing to population 2 after split.
    nu1: Final size of population 1 after growth.
    nu2: Size of population 2 after growth.
    Ts: The time between the split (=start of exp growth for pop1) and exp groth of pop1
    Tg1: time since the beggining of pop1 exp growth
    Tg2: Tg1+Tg2= time since the beggining of pop2 exp growth
    """
    nu1, nu2, s, Ts, Tg1, Tg2= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, Ts, nu1=1-s, nu2=s)
    # function for exp growth in pop 2 starting just at the split T
    nu2_func = lambda t: s*(nu2/s)**(t/(Tg2+Tg1))
    # pop1 immidiately starts to exponentialy grow after split (Ts+Tg1 unit of time ago)
    phi = Integration.two_pops(phi, xx, Tg2, nu1=1-s, nu2=nu2_func)
    # Then, Tg1 units of time ago, pop1 starts growing exponentialy 
    nu1_func = lambda t: (1-s)*(nu1/(1-s))**(t/Tg1)
    phi = Integration.two_pops(phi, xx, Tg1, nu1=nu1_func, nu2=nu2_func)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def split_bot_growth_rur_then_late_growth_urb_sym_mig(params, ns, pts):
    """
    Split into two populations, no migration. 
    Pop2 experiencing exponential growth few moments after split (from size 1-s).
    Pop2 starts its exponential growth later.
    s: proportion of Nref contributing to population 2 after split.
    nu1: Final size of population 1 after growth.
    nu2: Size of population 2 after growth.
    Ts: The time between the split (=start of exp growth for pop1) and exp groth of pop1
    Tg1: time since the beggining of pop1 exp growth
    Tg2: Tg1+Tg2= time since the beggining of pop2 exp growth
    """
    nu1, nu2, s, Ts, Tg1, Tg2, m= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, Ts, nu1=1-s, nu2=s, m12=m, m21=m)
    # function for exp growth in pop 2 starting just at the split T
    nu2_func = lambda t: s*(nu2/s)**(t/(Tg2+Tg1))
    # pop1 immidiately starts to exponentialy grow after split (Ts+Tg1 unit of time ago)
    phi = Integration.two_pops(phi, xx, Tg2, nu1=1-s, nu2=nu2_func, m12=m, m21=m)
    # Then, Tg1 units of time ago, pop1 starts growing exponentialy 
    nu1_func = lambda t: (1-s)*(nu1/(1-s))**(t/Tg1)
    phi = Integration.two_pops(phi, xx, Tg1, nu1=nu1_func, nu2=nu2_func, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def split_bot_growth_rur_then_late_growth_urb_asym_mig(params, ns, pts):
    """
    Split into two populations, no migration. 
    Pop2 experiencing exponential growth few moments after split (from size 1-s).
    Pop2 starts its exponential growth later.
    s: proportion of Nref contributing to population 2 after split.
    nu1: Final size of population 1 after growth.
    nu2: Size of population 2 after growth.
    Ts: The time between the split (=start of exp growth for pop1) and exp groth of pop1
    Tg1: time since the beggining of pop1 exp growth
    Tg2: Tg1+Tg2= time since the beggining of pop2 exp growth
    """
    nu1, nu2, s, Ts, Tg1, Tg2, m12, m21= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, Ts, nu1=1-s, nu2=s, m12=m12, m21=m21)
    # function for exp growth in pop 2 starting just at the split T
    nu2_func = lambda t: s*(nu2/s)**(t/(Tg2+Tg1))
    # pop1 immidiately starts to exponentialy grow after split (Ts+Tg1 unit of time ago)
    phi = Integration.two_pops(phi, xx, Tg2, nu1=1-s, nu2=nu2_func, m12=m12, m21=m21)
    # Then, Tg1 units of time ago, pop1 starts growing exponentialy 
    nu1_func = lambda t: (1-s)*(nu1/(1-s))**(t/Tg1)
    phi = Integration.two_pops(phi, xx, Tg1, nu1=nu1_func, nu2=nu2_func, m12=m12, m21=m21)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

################################################################

### Split, cst size then exp growth in pop2 the in pop1 ########

def split_bot_growth_urb_then_late_growth_rur_no_mig(params, ns, pts):
    """
    Split into two populations, no migration. 
    Pop2 experiencing exponential growth few moments after split (from size s).
    Pop1 starts its exponential growth later.
    s: proportion of Nref contributing to population 2 after split.
    nu1: Final size of population 1 after growth.
    nu2: Size of population 2 after growth.
    Ts: The time between the split (=start of exp growth for pop1) and exp groth of pop1
    Tg1: Tg1+Tg2= time since the beggining of pop1 exp growth
    Tg2: time since the beggining of pop2 exp growth
    """
    nu1, nu2, s, Ts, Tg1, Tg2= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, Ts, nu1=1-s, nu2=s)
    # function for exp growth in pop 2 starting at Tg1+Tg2
    nu1_func = lambda t: (1-s)*(nu1/(1-s))**(t/Tg1+Tg2)
    # pop1 immidiately starts to exponentialy grow after split
    phi = Integration.two_pops(phi, xx, Tg2, nu1=nu1_func, nu2=s)
    # Then, Tg2 units of time ago, pop1 starts growing exponentialy 
    nu2_func = lambda t: s*(nu2/s)**(t/(Tg2))
    phi = Integration.two_pops(phi, xx, Tg1, nu1=nu1_func, nu2=nu2_func)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def split_bot_growth_urb_then_late_growth_rur_sym_mig(params, ns, pts):
    """
    Split into two populations, no migration. 
    Pop2 experiencing exponential growth few moments after split (from size s).
    Pop1 starts its exponential growth later.
    s: proportion of Nref contributing to population 2 after split.
    nu1: Final size of population 1 after growth.
    nu2: Size of population 2 after growth.
    Ts: The time between the split (=start of exp growth for pop1) and exp groth of pop1
    Tg1: Tg1+Tg2= time since the beggining of pop1 exp growth
    Tg2: time since the beggining of pop2 exp growth
    """
    nu1, nu2, s, Ts, Tg1, Tg2, m= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, Ts, nu1=1-s, nu2=s, m12=m, m21=m)
    # function for exp growth in pop 2 starting at Tg1+Tg2
    nu1_func = lambda t: (1-s)*(nu1/(1-s))**(t/Tg1+Tg2)
    # pop1 immidiately starts to exponentialy grow after split
    phi = Integration.two_pops(phi, xx, Tg2, nu1=nu1_func, nu2=s, m12=m, m21=m)
    # Then, Tg2 units of time ago, pop1 starts growing exponentialy 
    nu2_func = lambda t: s*(nu2/s)**(t/(Tg2))
    phi = Integration.two_pops(phi, xx, Tg1, nu1=nu1_func, nu2=nu2_func, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def split_bot_growth_urb_then_late_growth_rur_asym_mig(params, ns, pts):
    """
    Split into two populations, no migration. 
    Pop2 experiencing exponential growth few moments after split (from size s).
    Pop1 starts its exponential growth later.
    s: proportion of Nref contributing to population 2 after split.
    nu1: Final size of population 1 after growth.
    nu2: Size of population 2 after growth.
    Ts: The time between the split (=start of exp growth for pop1) and exp groth of pop1
    Tg1: Tg1+Tg2= time since the beggining of pop1 exp growth
    Tg2: time since the beggining of pop2 exp growth
    """
    nu1, nu2, s, Ts, Tg1, Tg2, m12, m21= params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, Ts, nu1=1-s, nu2=s, m12=m12, m21=m21)
    # function for exp growth in pop 2 starting at Tg1+Tg2
    nu1_func = lambda t: (1-s)*(nu1/(1-s))**(t/Tg1+Tg2)
    # pop1 immidiately starts to exponentialy grow after split
    phi = Integration.two_pops(phi, xx, Tg2, nu1=nu1_func, nu2=s, m12=m12, m21=m21)
    # Then, Tg2 units of time ago, pop1 starts growing exponentialy 
    nu2_func = lambda t: s*(nu2/s)**(t/(Tg2))
    phi = Integration.two_pops(phi, xx, Tg1, nu1=nu1_func, nu2=nu2_func, m12=m12, m21=m21)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs
