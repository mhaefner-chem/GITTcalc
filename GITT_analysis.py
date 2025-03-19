# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 15:51:35 2025

@author: Anwender
"""

'''
function reads in the GITT data from a given file
TODO for OriginLab should read in data from sheet instead
'''
def fetch_GITT_data(file):
    data = {}
    print(file)
    
    labels = []
    
    splitter = ' '
    with open(file,mode='r') as f:
        lines = f.readlines()
        if ',' in lines[0]:
            splitter = ','
        elif '\t' in lines[0]:
            splitter = '\t'
        else:
            splitter = ' '
            
        required = ['time','Ewe','Capacity','spec. cap.']
            
        for line_index, line in enumerate(lines):
            if line_index == 0:
                for item in line.split(splitter):
                    item = item.split('\n')[0]
                    label = item
                    if item in ['Test Time (s)','time/s','Time (s)']:
                        label = 'time'
                    elif item in ['Voltage (V)','Ewe/V']:
                        label = 'Ewe'
                    elif item in ['Capacity (mAh)','Capacity/mA.h','Capacity/mAh']:
                        label = 'Capacity'
                    elif item in ['Special Capacity (mAh/g)']:
                        label = 'spec. cap.'
                        
                    unit = 1
                        
                    
                    data[label] = {"unit":unit,"values":[]}
                    labels.append(label)
            else:
                
                for item_index, item in enumerate(line.split(splitter)):
                    if labels[item_index] in required:
                        data[labels[item_index]]["values"].append(float(item))
                    
    return data

'''
Simple class to contain all settings
'''
class GITT_settings:
    def __init__(self, file):
        self.file_settings = file
        self.file_data ='test.dat'
        self.ion = 'Li'
        self.Z = 1
        self.m_AM = 1.118E-2    # g
        self.M = 97.871         # g/mol
        self.rho = 4.8          # g/cm³
        self.I0 = 1.49E-4       # A
        self.A_cont = 1.3273    # cm²
        self.theo_max_cap = 250 # mAh/g, spec. capacity
        self.c0_ion = 1         # starting stoichiometric Li/Na content of structure formula
        self.scale = 1.5
        self.limiter = 0.03
            
    
    '''
    Function to read settings for the system, measurement, and data analysis from file
    '''
    def fetch_GITT_settings(self):
        with open(self.file_settings, mode='r') as f:
            for line in f.readlines():
                if line.split()[0] == 'Z':
                    self.Z = int(line.split()[1])
                elif line.split()[0] == 'file':
                    self.file_data = line.split()[1]  
                elif line.split()[0] == 'theo_max_cap':
                    self.theo_max_cap = float(line.split()[1])
                elif line.split()[0] == 'm_AM':
                    self.m_AM = float(line.split()[1])
                elif line.split()[0] == 'M':
                    self.M = float(line.split()[1])
                elif line.split()[0] == 'rho':
                    self.rho = float(line.split()[1])
                elif line.split()[0] == 'I0':
                    self.I0 = float(line.split()[1])
                elif line.split()[0] == 'A_cont':
                    self.A_cont = float(line.split()[1])
                elif line.split()[0] == 'c0_ion':
                    self.c0_ion = float(line.split()[1])
                elif line.split()[0] == 'scale':
                    self.scale = float(line.split()[1])
                elif line.split()[0] == 'limiter':
                    self.limiter = float(line.split()[1])
                elif line.split()[0] == 'ion':
                    self.ion = line.split()[1]
        return
    
    




'''
function produces a numerical derivative of a given pair of x and y-values
slope at x determined with formula f(x+delta)-f(x-delta)/(2*delta)
'''
def get_numerical_derivative(x,y):
    
    derivative = []
    
    for i, value in enumerate(x):
        if i == 0:
            m = (y[i+1]-y[i])/(x[i+1]-x[i])
        elif i == len(x)-1:
            m = (y[i]-y[i-1])/(x[i]-x[i-1])
        else:
            m = (y[i+1]-y[i-1])/(x[i+1]-x[i-1])
        derivative.append(m)
        
    return derivative


if __name__ == '__main__':
    import numpy as np
    import math
    from scipy.stats import linregress
    
    has_cap = False
    has_specap = False
    
    # necessary values and constants for calculation for D and x_ion
    F = 96485.3321 # C/mol
    
    
    # retrieve data from given settings file
    # filename = 'settings_MFX_Arbin-P3.dat'
    # filename = '../settings_LNTO_AGR.dat'
    filename = '../settings_MFX.dat'
    settings = GITT_settings(filename)
    settings.fetch_GITT_settings()
                                            
    GITT_data = fetch_GITT_data(settings.file_data)
    
    x = GITT_data['time']['values']
    y = GITT_data['Ewe']['values']
    
    x_spec_cap = []
    x_ion = []
    y_cycle = []
    
    if 'Capacity' in GITT_data:
        x_cap = GITT_data['Capacity']['values']
        has_cap = True
    elif 'spec. cap.' in GITT_data:
        x_spec_cap = GITT_data['spec. cap.']['values']
        for value in x_spec_cap:
            x_ion.append(-value/settings.theo_max_cap + settings.c0_ion)
            y_cycle.append(0)
        has_specap = True
    else:
        x_capacity = [0]*len(x)
    y_deriv = get_numerical_derivative(x, y)
    
    
    # calculate specific capacity and ion content at every given time
    charge = True
    current_cycle = 0 # charge-discharge cycle number
    ref_cap = 0
    if has_cap:
        y_cycle = []
        # real_max_cap = max(x_cap)
        for i, cap in enumerate(x_cap):    
            y_cycle.append(current_cycle)
            if i > 1 and x_cap[i]-x_cap[i-1] < 0:
                current_cycle += 1
                if charge:
                    charge = False
                    ref_cap = ref_cap+x_cap[i-1]
                else:
                    charge = True
                    ref_cap = ref_cap-x_cap[i-1]
                
            
        
            if charge == True:
                spec_cap = (ref_cap+cap)/settings.m_AM
            else:
                spec_cap = (ref_cap-cap)/settings.m_AM
            
            x_spec_cap.append(spec_cap)
            x_ion.append(-spec_cap/settings.theo_max_cap + settings.c0_ion)
    
    
    
    # get special points
    current_on = []
    current_off = []
    on_times = []
    off_times = []
    system_equil = []
    
    # determine the sensitivity of detecting the derivative jumps
    
    # scales the check to the mean of the absolute values of the derivative
    y_deriv_cutoff = settings.limiter*np.mean(list(map(abs, y_deriv)))

    load = False  
    for i, value in enumerate(y_deriv):
        # detects positive derivative jump
        if y_deriv[i] > abs(settings.scale*y_deriv[i-1]) and y_deriv[i] > y_deriv_cutoff and load == False:
            # switches meaning of jump depending of whether the cell is currently being charged or discharged
            if len(current_on) == 0 or y[current_on[-1]] < y[i]:
                current_on.append(i)
                on_times.append(x[i])
            else:
               current_off.append(i)
               off_times.append(x[i])
            load = True
        # detects negative derivative jump, same logic as for positive derivative jump but flipped
        elif y_deriv[i] < -abs(settings.scale*y_deriv[i-1]) and y_deriv[i] < -y_deriv_cutoff and load == True:
            if len(current_on) == 0 or y[current_on[-1]] < y[i]:
                current_off.append(i)
                off_times.append(x[i])
            else:
               current_on.append(i)
               on_times.append(x[i])
            load = False
    
    # evaluate E1-E4, tau
    GITT_results = []
    D_Li = []
    D_x = []
    D_cycle = []
    D_spec_cap = []
    D_voltage = []
    
    off = 0
    for i, on in enumerate(current_on):
        # if i > 1 and y_cycle[current_on[i+1]] > 1:
            # break
        
        for off in current_off:
            if off > on:
                break
        # for equil in system_equil:
        #     if equil > on:
        #         break
        if i == len(current_on)-1:
            break
        
        
        E1  = y[on]
        E3  = y[off]
        E4  = y[current_on[i+1]]
        tau = x[off]-x[on]
        # print(x[on],x[off],y[on],y[off])
        
        # linear regression for sqrt-relationship for E2
        interval_x = []
        interval_y = []
        for j in range(on+int((on-off)/2),off):
            interval_x.append(np.sqrt(x[j]))
            interval_y.append(y[j])
        try:
            regress_param = linregress(interval_x,interval_y)
        except:
            continue
        m = [regress_param.slope,regress_param.stderr]
        b = [regress_param.intercept,regress_param.intercept_stderr]
        E2 = m[0]*np.sqrt(x[on]) + b[0]
        
        GITT_results.append((E1,E2,E3,E4,tau,x[on],x[on],x[off],regress_param.rvalue**2))
        if has_cap or has_specap:
            D_Li.append(x_ion[current_on[i]])
            D_spec_cap.append(x_spec_cap[current_on[i]])
            D_cycle.append(y_cycle[current_on[i]])
            
        D_x.append(x[current_on[i]])
        D_voltage.append(y[current_on[i]])
        

    V_mol = settings.M/settings.rho
    print(V_mol)
    
    D = []

    for result in GITT_results:
        D.append(4/(np.pi*result[4]) * (settings.m_AM * V_mol/(settings.M*settings.A_cont))**2 * ((result[3]-result[0])/(result[2]-result[1]))**2)

        # print(result,D)
        # print(result[4])
        
    with open('output.csv',mode='w') as f:
        if has_cap or has_specap:
            f.write('{:16},{:16},{:16},{:16},{:16},{:16}\n'.format('time/s','Ewe/V','x_{}'.format(settings.ion),'SpecCap/mAh/g','D/cm^2/s','Cycle'))
            for i,value in enumerate(D):
                f.write('{:16.10e},{:16.10e},{:16.10e},{:16.10e},{:16.10e},{:4}\n'.format(D_x[i],D_voltage[i],D_Li[i],D_spec_cap[i],D[i],D_cycle[i]))
        else:
            f.write('{:16},{:16},{:16}\n'.format('time/s','Ewe/V','D/cm^2/s'))
            for i,value in enumerate(D):
                f.write('{:16.10e},{:16.10e},{:16.10e}\n'.format(D_x[i],D_voltage[i],D[i]))
    
    
    # everything after here is only for plotting the results
    plot = True
    if plot == True:
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        fs = 16
        fig, ax = plt.subplots()
        colors = ["green","blue","red"]
        alphas = [1.0,1.0,1.0]
        labels = ['E1','E2','E3']
        for i in range(3):
            tmp = [[],[]]
            for result in GITT_results:
                tmp[0].append(result[i+5])
                tmp[1].append(result[i])
            ax.scatter(tmp[0],tmp[1],marker="x",color=colors[i],zorder=50,label=labels[i],alpha=alphas[i])
        
        ax.plot(x,y,linestyle='-',label='E')
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
            
        
        # tmpx = np.linspace(62000,63000,100)
        # tmpy = []
        # for value in tmpx:
        #     tmpy.append(m[0]*np.sqrt(value)+b[0])
        # ax.plot(tmpx,tmpy,color="red",linestyle="-")
        
        ax2 = ax.twinx()
        # ax2.plot(x,y_deriv, color = 'black',linewidth=1, zorder=0,linestyle="-",label="dE/dt")  
        ax2.scatter(D_x,D,color="red",marker="x",label="D")
        ax2.set_ylim(1E-11,5E-9)
        ax2.axhline(0,color="#000000",linewidth=1)
        
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax2.legend(h1+h2, l1+l2,loc=1,fontsize=fs)
        
        # ax.set_xlim(610000,640000)
        # ax.set_ylim(4.17,4.32)
        ax2.set_yscale('log')
        # ax2.set_ylim(-0.001,0.001)
        ax2.set_ylim(0,1.5E-9)
        plt.yticks(fontsize=fs)
        
        ax.grid(zorder=-50,linestyle="--",alpha=0.66)
        ax.set_xlabel('Time t ($10^6$ s)',fontsize=fs) 
        ax.set_ylabel('Voltage E (V)',fontsize=fs) 
        # ax2.set_ylabel('Derivative dE/dt [V/s]',fontsize=fs) 
        ax2.set_ylabel('Diffusion rate D ($10^{-9}$ cm²/s)',fontsize=fs) 
        ax.set_title("GITT Plot\nDiffusion and Voltage against Time",fontsize=fs+4)
        
        plt.show()
        # fig.savefig("GITT.png", format='png',bbox_inches="tight",dpi=300)

        if has_cap or has_specap:
            fig, ax = plt.subplots()
    
            # cmap = mpl.cm.viridis
            cmap = mpl.cm.tab10
            values = {}
            for i,value in enumerate(D):
                
                if D_cycle[i]%2 == 0:
                    label = "charge cycle "
                else:
                    label = "discharge cycle "
                    
                label += str(int(D_cycle[i]/2)+1)
                    
                if i == 0:
                    values[D_cycle[i]]=[[],[],label]
                elif D_cycle[i] != D_cycle[i-1]:
                    values[D_cycle[i]]=[[],[],label]
    
                values[D_cycle[i]][0].append(D_Li[i])
                values[D_cycle[i]][1].append(D[i])
    
            for key,item in values.items():
                ax.scatter(item[0],item[1],marker='x',label=item[2],color=cmap(key/10),s=25,zorder=50)
            
            # make charge filled, discharge empty markers!!!
            
            # ax.scatter(x_ion,y,marker='o',label='E',c="black",s=0.1,zorder=50)
            # ax.scatter(D_Li,D,marker='x',label='D',c=color,cmap=cmap,s=25,zorder=50)
            # ax.scatter(x,x_ion,marker='x',label='E',c="red",s=1)
            # cbar = fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap),
                         # ax=ax, orientation='vertical', label='')
            ax.grid(zorder=-50,linestyle="--",alpha=0.66)
            ax.axhline(0,color="#000000",linewidth=1)
            # ax.set_ylim(0,0.2E-8)
            ax.set_ylim(1E-11,1E-8)
            ax.set_yscale('log')
            ax.legend(fontsize=fs)
            plt.xticks(fontsize=fs)
            plt.yticks(fontsize=fs)
            ax.set_xlabel('x$_{\mathregular{Na}}$',fontsize=fs) 
            ax.set_ylabel('D ($10^{-9}$ cm²/s)',fontsize=fs) 
            ax.set_title("GITT plot\n Diffusion at Varying x$_{\mathregular{Na}}$\n(x$_{\mathregular{Na}}$ at beginning of each titration step)",fontsize=fs)
            plt.show()
        
    
    
    
    
    