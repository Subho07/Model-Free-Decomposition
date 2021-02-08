from osgeo import gdal
import numpy as np
import warnings
warnings.filterwarnings("ignore")

def conv2d(a, f):
    filt = np.zeros(a.shape)
    wspad = int(f.shape[0]/2)
    s = f.shape + tuple(np.subtract(a.shape, f.shape) + 1)
    strd = np.lib.stride_tricks.as_strided
    subM = strd(a, shape = s, strides = a.strides * 2)
    filt_data = np.einsum('ij,ijkl->kl', f, subM)
    filt[wspad:wspad+filt_data.shape[0],wspad:wspad+filt_data.shape[1]] = filt_data
    return filt

def read_bin(file):
    ds = gdal.Open(file)
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray()
   
    return arr

def write_bin(file,wdata,refData):
                
    ds = gdal.Open(refData)
    [cols, rows] = wdata.shape
            
    driver = gdal.GetDriverByName("ENVI")
    outdata = driver.Create(file, rows, cols, 1, gdal.GDT_Float32)
    outdata.SetGeoTransform(ds.GetGeoTransform())##sets same geotransform as input
    outdata.SetProjection(ds.GetProjection())##sets same projection as input
                
    outdata.SetDescription(file)
    outdata.GetRasterBand(1).WriteArray(wdata)
                
    outdata.FlushCache() ##saves to disk!! 
    
def mf3cc_powers(C3_folder,ws, chi_in = -45):

    C11 = read_bin(C3_folder+"/C11.bin")
    C22 = read_bin(C3_folder+"/C22.bin")

    C12_i = read_bin(C3_folder+'/C12_imag.bin')
    C12_r = read_bin(C3_folder+'/C12_real.bin')

    C12 = C12_r + 1j*C12_i
        
    c11_T1 = C11
    c12_T1 = C12
    c21_T1 = np.conj(C12)
    c22_T1 = C22
                
                
    kernel = np.ones((ws,ws),np.float32)/(ws*ws)

            
    c11_T1r = conv2d(np.real(c11_T1),kernel)
    c11_T1i = conv2d(np.imag(c11_T1),kernel)
    c11s = c11_T1r+1j*c11_T1i

    c12_T1r = conv2d(np.real(c12_T1),kernel)
    c12_T1i = conv2d(np.imag(c12_T1),kernel)
    c12s = c12_T1r+1j*c12_T1i

    c21_T1r = conv2d(np.real(c21_T1),kernel)
    c21_T1i = conv2d(np.imag(c21_T1),kernel)
    c21s = c21_T1r+1j*c21_T1i


    c22_T1r = conv2d(np.real(c22_T1),kernel)
    c22_T1i = conv2d(np.imag(c22_T1),kernel)
    c22s = c22_T1r+1j*c22_T1i

    c2_det = (c11s*c22s-c12s*c21s)
    c2_trace = c11s+c22s
    m1 = np.real(np.sqrt(1.0-(4.0*c2_det/np.power(c2_trace,2))))

    s0 = c11s + c22s
    s1 = c11s - c22s
    s2 = (c12s + c21s)

    if (chi_in >= 0):
        s3 = (1j*(c12s - c21s))
    if (chi_in < 0):
        s3 = -(1j*(c12s - c21s))
                
    SC = ((s0)-(s3))/2
    OC = ((s0)+(s3))/2

    h = (OC-SC)
    span = c11s + c22s

    val = ((m1*s0*h))/((SC*OC + (m1**2)*(s0**2)))
    thet = np.real(np.arctan(val))
    theta_CP = np.rad2deg(thet)

    Ps_CP= (((m1*(span)*(1.0+np.sin(2*thet))/2)))
    Pd_CP= (((m1*(span)*(1.0-np.sin(2*thet))/2)))
    Pv_CP= (span*(1.0-m1))


    infile = C3_folder+'/C11.bin'

    ofilegrvi = C3_folder+'/Theta_CP.bin'
    write_bin(ofilegrvi,theta_CP,infile)
                
    ofilegrvi1 = C3_folder+'/Pd_CP.bin'
    write_bin(ofilegrvi1,Pd_CP,infile)
                
    ofilegrvi2 = C3_folder+'/Ps_CP.bin'
    write_bin(ofilegrvi2,Ps_CP,infile)
                
    ofilegrvi3 = C3_folder+'/Pv_CP.bin'
    write_bin(ofilegrvi3,Pv_CP,infile)

    return Ps_CP, Pd_CP, Pv_CP
