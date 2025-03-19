# -*- coding: utf-8 -*-
'''
Created on Sat Feb  8 15:51:35 2025

@author: mhaefner-chem
'''

'''
function reads in the GITT data from a given file
INFO: this needs to be linked to OriginLab
Format: 
    dictionary 'data', keys from list ['time','volt','cap','spec_cap']
    time (continuous measurement time) and volt (measured voltage) are required
    either cap (capacity) or spec_cap (special capacity) are requied
    
    datapoints are stored in lists:
    data[key] = []
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
            
        required = ['time','volt','cap','spec_cap']
            
        for line_index, line in enumerate(lines):
            if line_index == 0:
                for item in line.split(splitter):
                    item = item.split('\n')[0]
                    label = item
                    if item in ['Test Time (s)','time/s','Time (s)']:
                        label = 'time'
                    elif item in ['Voltage (V)','Ewe/V']:
                        label = 'volt'
                    elif item in ['Capacity (mAh)','Capacity/mA.h','Capacity/mAh']:
                        label = 'cap'
                    elif item in ['Special Capacity (mAh/g)']:
                        label = 'spec_cap'
                        
                    data[label] = []
                    labels.append(label)
            else:
                
                for item_index, item in enumerate(line.split(splitter)):
                    if labels[item_index] in required:
                        data[labels[item_index]].append(float(item))
                    
    return data

'''
Simple class to contain all settings
INFO: this needs to be linked to OriginLab
Format:
    settings are stored as variables in the class
'''
class GITT_settings:
    '''
    Function to initialize settings for the system, measurement, and data analysis
    '''
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
    
'''
function outputs diffusion data
INFO: this needs to be linked to OriginLab
Format: 
    dictionary 'data_out', keys from list ['time','volt','spec_cap','cycle','ion','diff']
    'time' is measurement time and 'volt' is measured voltage associated with diffusion rate 'diff'
    
    if capacity data is supplied, output also contains
    'spec_cap' (special capacity)
    'ion' (content of diffusing ion, e.g., Li or Na)
    'cycle' (current cycle, usually 0 is first charge, 1 is first discharge, 2 is second charge, etc.)
    
    datapoints are stored in lists:
    data_out[key] = []
'''
def write_GITT_data(data_out,has_capacity):
    
    with open('../output.csv',mode='w') as f:
        if has_capacity:
            f.write('{:16},{:16},{:16},{:16},{:16},{:16}\n'.format('time/s','volt/V','x_{}'.format(settings.ion),'SpecCap/mAh/g','D/cm^2/s','Cycle'))
            for i,value in enumerate(data_out['diff']):
                f.write('{:16.10e},{:16.10e},{:16.10e},{:16.10e},{:16.10e},{:4}\n'.format(data_out['time'][i],data_out['volt'][i],data_out['ion'][i],data_out['spec_cap'][i],data_out['diff'][i],data_out['cycle'][i]))
        else:
            f.write('{:16},{:16},{:16}\n'.format('time/s','volt/V','D/cm^2/s'))
            for i,value in enumerate(data_out['diff']):
                f.write('{:16.10e},{:16.10e},{:16.10e}\n'.format(data_out['time'][i],data_out['volt'][i],data_out['diff'][i]))
                
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

'''
function for plotting results from analyzed GITT data
messy code, refactoring not yet planned
'''
def plot_results(GITT_refined,D_out):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    fs = 16
    fig, ax = plt.subplots()
    colors = ["green","blue","red"]
    alphas = [1.0,1.0,1.0]
    labels = ['E1','E2','E3']
    for i in range(3):
        tmp = [[],[]]
        for result in GITT_refined:
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
    ax2.scatter(D_out['time'],D_out['diff'],color="red",marker="x",label="D")
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
        for i,value in enumerate(D_out['diff']):
            
            if D_out['cycle'][i]%2 == 0:
                label = "charge cycle "
            else:
                label = "discharge cycle "
                
            label += str(int(D_out['cycle'][i]/2)+1)
                
            if i == 0:
                values[D_out['cycle'][i]]=[[],[],label]
            elif D_out['cycle'][i] != D_out['cycle'][i-1]:
                values[D_out['cycle'][i]]=[[],[],label]

            values[D_out['cycle'][i]][0].append(D_out['ion'][i])
            values[D_out['cycle'][i]][1].append(D_out['diff'][i])

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

if __name__ == '__main__':
    import numpy as np
    from scipy.stats import linregress
    
    # retrieve settings for calculation and file containing GITT data
    # INFO: should not be necessary for OriginLab implementation
    filename = '../settings_MFX_arbin.dat'
    settings = GITT_settings(filename)
    settings.fetch_GITT_settings()
                      
    # retrieve GITT data, naturally works with Biologic output format                
    GITT_data = fetch_GITT_data(settings.file_data)
    
    # initial data transformation, time as x-axis, voltage as y-axis
    x = GITT_data['time']
    y = GITT_data['volt']
    
    # alternative x-axis for special capacity and ion content
    # alternative y-axis for current cycle

    
    # assigns capacity if it was provided
    # prefers special capacity over pure capacity, if both are provided
    has_cap = False
    has_specap = False
    if 'spec_cap' in GITT_data:
        x_spec_cap = GITT_data['spec_cap']
        has_specap = True
    elif 'cap' in GITT_data:
        x_cap = GITT_data['cap']
        has_cap = True
    else:
        x_capacity = [0]*len(x)    
    
    # calculate specific capacity, current cycle, and ion content at every given time
    charge = True
    current_cycle = 0 # charge-discharge cycle number
    ref_cap = 0
    x_ion = []
    y_cycle = []
    # this logic only works with Biologic output!!!
    # i.e., capacity is reset to 0 for every new charge or discharge cycle
    # capacity only ever rises and never decreases
    if has_cap:
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
            
    # this logic works with capacity that rises during charge and decreases during discharge
    # i.e., capacity goes up for charge and down for discharges
    elif has_specap:
        charge = True
        for i, spec_cap in enumerate(x_spec_cap):
            x_ion.append(-spec_cap/settings.theo_max_cap + settings.c0_ion)
            y_cycle.append(current_cycle)
            if i > 1:
                if x_spec_cap[i]-x_spec_cap[i-1] < 0 and charge == True:
                    current_cycle += 1
                    charge = False
                elif x_spec_cap[i]-x_spec_cap[i-1] > 0 and charge == False:
                    current_cycle += 1
                    charge = True
    
    # get numerical derivative of voltage
    # cutoff determines minimum jump in derivative required for it to be counted
    y_deriv = get_numerical_derivative(x, y)
    y_deriv_cutoff = settings.limiter*np.mean(list(map(abs, y_deriv)))

    # this part detects when the current is applied and removed
    # makes this less dependent on format of GITT data
    current_on = []
    current_off = []
    on_times = []
    off_times = []
    system_equil = []
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
    
    # evaluate E1-E4, charging time tau
    D_out = {
        'ion':      [],
        'spec_cap': [],
        'cycle':    [],
        'time':     [],
        'volt':     [],
        'diff':     []
        }
    GITT_refined = [] # ONLY REQUIRED FOR PLOTTING
    
    V_mol = settings.M/settings.rho    
    F = 96485.3321 # C/mol
    off = 0
    for i, on in enumerate(current_on):    
        # determines the next point after on at which current is turned off
        for off in current_off:
            if off > on:
                break
        
        # determination of E1, E3, and tau
        E1 = y[on]
        E3 = y[off]
        tau = x[off]-x[on]
        
        # E4 requires logic to properly treat final titration
        if i < len(current_on)-1:
            E4 = y[current_on[i+1]]
        else:
            E4 = y[-1]
            
        # E2 requires linear regression for sqrt-behavior while current is applied
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
        
        GITT_refined.append((E1,E2,E3,E4,tau,x[on],x[on],x[off])) # ONLY REQUIRED FOR PLOTTING
        
        # collect data for output
        if has_cap or has_specap:
            D_out['ion'].append(x_ion[current_on[i]])
            D_out['spec_cap'].append(x_spec_cap[current_on[i]])
            D_out['cycle'].append(y_cycle[current_on[i]])   
        D_out['time'].append(x[current_on[i]])
        D_out['volt'].append(y[current_on[i]])
        
        # calculation of the diffusion constants
        dE_s = E4-E1
        dE_t = E3-E2
        D = 4/(np.pi*tau) * (settings.m_AM * V_mol/(settings.M*settings.A_cont))**2 * ((dE_s)/(dE_t))**2
        D_out['diff'].append(D) 
    
    # writes diffusion data derived from GITT measurement
    write_GITT_data(D_out,has_cap or has_specap)

    # for plotting purposes only
    plot = True
    if plot == True:
        plot_results(GITT_refined,D_out)
        