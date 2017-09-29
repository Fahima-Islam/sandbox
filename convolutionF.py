#matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
import math
def convolve(Iin,k,opt):
	#str_arr_I = raw_input('insert only the values of the first 1D array with the space between them:').split(' ') 
	#Iin=0.0* np.ones(len(str_arr_I))
	#for i,j in zip (str_arr_I, xrange(len(str_arr_I))):    
   		# Iin[j]=i
	#str_arr_k = raw_input('insert only values of the second 1D array with the space between them:').split(' ') 
	#k=0.0* np.ones(len(str_arr_k))
	#for i,j in zip (str_arr_k, xrange(len(str_arr_k))):    
   		# k[j]=i

	kT=k[::-1]   
	if len(Iin)<len(k):
    		kT=Iin[::-1]    
	kT0=kT
	kT1=kT
	if len(k)!=len(Iin): 
    		kT_=[0]*abs(len(k)-len(Iin))
    		kT0=np.concatenate((kT_,kT), axis=0)
    		kT1=np.concatenate((kT,kT_), axis=0)
    
	length_s=0
	if len(Iin)>len(k):
    		length_s=len(Iin)
	else:
    		length_s=len(k)
    
	if len(Iin)<len(k):
    		Iin=k
    
	M,T=np.meshgrid(Iin,Iin)
	il=np.tril(M[:(len(kT)-1),:])
	iu=np.triu(M)
	mask = il>0
	mask_1 = iu>0
	justified_mask = np.sort(mask,1)
	justified_mask_1 = np.sort(mask_1,1)
	justified_mask = justified_mask[:,::]
	justified_mask_1 = justified_mask_1[:,::-1]
	out = np.zeros_like(il[:,:]) 
	out_1 = np.zeros_like(iu)
	out[justified_mask] = il[:,:][mask]
	out_1[justified_mask_1] = iu[mask_1]
	Rr=np.concatenate((np.dot(out,kT0),np.dot(kT1,out_1)),axis=0)
	#print("for full mode {}".format(Rr) )
	length_f=len(Rr)
	off=length_f-length_s
	off_eachend=off/2
	f_idx=int(off_eachend)
	e_idx=f_idx+length_s
	Sr=Rr[f_idx:e_idx]
	#print("for same mode {}".format(Sr))
	if opt=='full':	
		return(Rr)
	if opt=='same':
		return(Sr)

import convolutionF as F
def deconvolve(sig,mask,deconV,conv,option,value):
	sig0=sig
	mask_mir=mask[::-1]
	m_tst=F.convolve(sig,mask,conv)
	deconv = deconV
	
	def main(deconv,mask,sig0,mask_mir,conv):
		sigC=F.convolve(deconv,mask,conv)
		relative_blur=sig0/sigC
		with np.errstate(divide='ignore'):
			relative_blur[np.isinf(relative_blur)] = -2
		deconvP=deconv*F.convolve(relative_blur,mask_mir,conv)
		error=np.abs(deconvP-deconv)
		deconv=deconvP
		return(deconv,error)

	if option=='iteration':
		error=0		
		for i in xrange(value):
			deconv,error=main(deconv,mask,sig0,mask_mir,conv)
    			#sigC=F.convolve(deconv,mask,conv)
    			#relative_blur=sig0/sigC
    			#deconvP=deconv*F.convolve(relative_blur,mask_mir,conv)
			#error=np.abs(deconvP-deconv)
                        #deconv=deconvP
		#print('error achieved: {}'.format(error))

	if option=='error':
		it=0
		while True:
			deconv,error=main(deconv,mask,sig0,mask_mir,conv)
    			#sigC=F.convolve(deconv,mask,conv)
    			#relative_blur=sig0/sigC
    			#deconvP=deconv*F.convolve(relative_blur,mask_mir,conv)
    			#error=np.abs(deconvP-deconv)
    			#deconv=deconvP
			it=it+1
			#print('number of iteration: {}'.format(it))
			if np.all(error<value):
				break
		print('number of iteration: {}'.format(it))
	return(deconv)

import convolutionF as F
def deconvolve_TV(sig,mask,deconV,conv,eps,rgP,option,value):
        sig0=sig
        mask_mir=mask[::-1]
        m_tst=F.convolve(sig,mask,conv)
        deconv = deconV

        def main(deconv,mask,sig0,mask_mir,conv,eps,rgP):
                sigC=F.convolve(deconv,mask,conv)
                relative_blur=sig0/sigC
                with np.errstate(divide='ignore'):
                        relative_blur[np.isinf(relative_blur)] = -2
		grad=np.gradient(deconv)
		#grad=np.diff(deconv)
		#norm=np.linalg.norm(grad, ord=1)
		norm=np.sqrt(grad**2)
		mod_norm=np.sqrt(eps**2+norm**2)
		division=(grad)/mod_norm
		division[np.isnan(division)] = 0.0
		with np.errstate(divide='ignore'):
			division[np.isinf(division)] = -2
		#divergence=np.sum(np.gradient(division))
		divergence=np.gradient(division)
		div_rgp=rgP*divergence
		#if np.any(div_rgp>1):
		#mask=(div_rgp)>1
		#div_rgp[:][mask]=divergence[mask]
		#mask2=(1-div_rgp)<0.02
		#div_rgp[:][mask2]=(rgP/3)*(divergence[mask2])
		#div_rgp[np.any(abs(div_rgp)>1)]=divergence
			#rgp=divergence/divergence
		#div_rgp=rgP*division
		deconvP=(deconv/(1-(div_rgp)))*F.convolve(relative_blur,mask_mir,conv)
		error=np.abs(deconvP-deconv)
		deconv=deconvP
		return(deconv,error)

        if option=='iteration':
                error=0         
                for i in xrange(value):
                        deconv,error=main(deconv,mask,sig0,mask_mir,conv,eps,rgP)
                        #sigC=F.convolve(deconv,mask,conv)
                        #relative_blur=sig0/sigC
                        #deconvP=deconv*F.convolve(relative_blur,mask_mir,conv)
                        #error=np.abs(deconvP-deconv)
                        #deconv=deconvP
                #print('error achieved: {}'.format(error))

        if option=='error':
                it=0
                while True:
                        deconv,error=main(deconv,mask,sig0,mask_mir,conv,eps,rgP)
                        #sigC=F.convolve(deconv,mask,conv)
                        #relative_blur=sig0/sigC
                        #deconvP=deconv*F.convolve(relative_blur,mask_mir,conv)
                        #error=np.abs(deconvP-deconv)
                        #deconv=deconvP
                        it=it+1
                        #print('number of iteration: {}'.format(it))
                        if np.all(error<value):
                                break
                print('number of iteration: {}'.format(it))
        return(deconv)

