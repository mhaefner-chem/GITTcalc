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
                    print(label)
                    if item in ['Test Time (s)','time/s','Time (s)']:
                        label = 'time'
                    elif item in ['Voltage (V)','Ewe/V']:
                        label = 'volt'
                    elif item in ['Capacity (mAh)','Capacity/mA.h','Capacity/mAh']:
                        label = 'cap'
                    elif item in ['Specific Capacity (mAh/g)']:
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
# class GITT_settings:
#     '''
#     Function to initialize settings for the system, measurement, and data analysis
#     '''
#     def __init__(self):
        
#         self.file_settings = ''
#         self.ion = 'Li'
#         self.Z = 1
#         self.m_AM = 1.118E-2    # g
#         self.M = 97.871         # g/mol
#         self.rho = 4.8          # g/cm³
#         self.I0 = 1.49E-4       # A
#         self.A_cont = 1.3273    # cm²
#         self.theo_max_cap = 250 # mAh/g, spec. capacity
#         self.c0_ion = 1         # starting stoichiometric Li/Na content of structure formula
#         self.scale = 1.5
#         self.limiter = 0.03
#         self.plot = {
#             'DV_t':True,
#             'D_x':False
#             }
#         self.has_cap = {
#             'cap':False,
#             'spec_cap':False
#             }
    
#     '''
#     Function to read settings for the system, measurement, and data analysis from file
#     '''
#     def fetch_GITT_settings(self):
        
#         with open(self.file_settings, mode='r') as f:
#             for line in f.readlines():
#                 if line.split()[0] == 'Z':
#                     self.Z = int(line.split()[1])
#                 elif line.split()[0] == 'file':
#                     self.file_data = line.split()[1]  
#                 elif line.split()[0] == 'theo_max_cap':
#                     self.theo_max_cap = float(line.split()[1])
#                 elif line.split()[0] == 'm_AM':
#                     self.m_AM = float(line.split()[1])
#                 elif line.split()[0] == 'M':
#                     self.M = float(line.split()[1])
#                 elif line.split()[0] == 'rho':
#                     self.rho = float(line.split()[1])
#                 elif line.split()[0] == 'I0':
#                     self.I0 = float(line.split()[1])
#                 elif line.split()[0] == 'A_cont':
#                     self.A_cont = float(line.split()[1])
#                 elif line.split()[0] == 'c0_ion':
#                     self.c0_ion = float(line.split()[1])
#                 elif line.split()[0] == 'scale':
#                     self.scale = float(line.split()[1])
#                 elif line.split()[0] == 'limiter':
#                     self.limiter = float(line.split()[1])
#                 elif line.split()[0] == 'ion':
#                     self.ion = line.split()[1]
    
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
def write_GITT_data(savefile,data_out,settings):
    
    ion = settings['ion'].entry.get()
    with savefile as f:
        if settings['cap'] or settings['spec_cap']:            
            f.write('{:16},{:16},{:16},{:16},{:16},{:16}\n'.format('time/s','volt/V','x_{}'.format(ion),'SpecCap/mAh/g','D/cm^2/s','Cycle'))
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

class plot_window:
    # initializes the window and default plotting data
    def __init__(self,GITT_data,GITT_refined,D_out,settings):
        
        self.root = create_window("1000x700+120+120", "Atomic Form Factor Plot")
            
        self.data = GITT_data
        self.refined = GITT_refined
        self.D = D_out
        self.settings = settings
        
        self.dpi_default = 100
        self.dpi_set = str(self.dpi_default)
        
        self.draw_window()
        
    # populates the window with widgets
    def draw_window(self):
        self.plot_form_factors()
    
    
    # creates the plot with matplotlib
    def plot_form_factors(self): 
        
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,  
        NavigationToolbar2Tk)
        
        fs = 16
        
        DV_t = True
        D_x = False
        
        if DV_t:
            print(123)
            fig, ax = plt.subplots()
            colors = ["green","blue","red"]
            alphas = [1.0,1.0,1.0]
            labels = ['E1','E2','E3']
            for i in range(3):
                tmp = [[],[]]
                for result in self.refined:
                    tmp[0].append(result[i+5])
                    tmp[1].append(result[i])
                ax.scatter(tmp[0],tmp[1],marker="x",color=colors[i],zorder=50,label=labels[i],alpha=alphas[i])
            
            ax.plot(self.data['time'],self.data['volt'],linestyle='-',label='E',marker='x')
            plt.xticks(fontsize=fs)
            plt.yticks(fontsize=fs)
                
            
            # tmpx = np.linspace(62000,63000,100)
            # tmpy = []
            # for value in tmpx:
            #     tmpy.append(m[0]*np.sqrt(value)+b[0])
            # ax.plot(tmpx,tmpy,color="red",linestyle="-")
            
            ax2 = ax.twinx()
            # ax2.plot(x,y_deriv, color = 'black',linewidth=1, zorder=0,linestyle="-",label="dE/dt")  
            ax2.scatter(self.D['time'],self.D['diff'],color="red",marker="x",label="D")
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

        # if settings['cap'] or settings['spec_cap'] and D_x:
        #     fig, ax = plt.subplots()

        #     # cmap = mpl.cm.viridis
        #     cmap = mpl.cm.tab10
        #     values = {}
        #     for i,value in enumerate(D_out['diff']):
                
        #         if D_out['cycle'][i]%2 == 0:
        #             label = "charge cycle "
        #         else:
        #             label = "discharge cycle "
                    
        #         label += str(int(D_out['cycle'][i]/2)+1)
                    
        #         if i == 0:
        #             values[D_out['cycle'][i]]=[[],[],label]
        #         elif D_out['cycle'][i] != D_out['cycle'][i-1]:
        #             values[D_out['cycle'][i]]=[[],[],label]

        #         values[D_out['cycle'][i]][0].append(D_out['ion'][i])
        #         values[D_out['cycle'][i]][1].append(D_out['diff'][i])

        #     for key,item in values.items():
        #         ax.scatter(item[0],item[1],marker='x',label=item[2],color=cmap(key/10),s=25,zorder=50)
            
        #     # make charge filled, discharge empty markers!!!
            
        #     # ax.scatter(x_ion,y,marker='o',label='E',c="black",s=0.1,zorder=50)
        #     # ax.scatter(D_Li,D,marker='x',label='D',c=color,cmap=cmap,s=25,zorder=50)
        #     # ax.scatter(x,x_ion,marker='x',label='E',c="red",s=1)
        #     # cbar = fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap),
        #                  # ax=ax, orientation='vertical', label='')
        #     ax.grid(zorder=-50,linestyle="--",alpha=0.66)
        #     ax.axhline(0,color="#000000",linewidth=1)
        #     # ax.set_ylim(0,0.2E-8)
        #     ax.set_ylim(1E-11,1E-8)
        #     ax.set_yscale('log')
        #     ax.legend(fontsize=fs)
        #     plt.xticks(fontsize=fs)
        #     plt.yticks(fontsize=fs)
        #     ax.set_xlabel('x$_{\mathregular{Na}}$',fontsize=fs) 
        #     ax.set_ylabel('D ($10^{-9}$ cm²/s)',fontsize=fs) 
        #     ax.set_title("GITT plot\n Diffusion at Varying x$_{\mathregular{Na}}$\n(x$_{\mathregular{Na}}$ at beginning of each titration step)",fontsize=fs)
        

            
        # creates and places Tkinter canvas for the matplotlib figure
        canvas = FigureCanvasTkAgg(fig, master = self.root)   
        canvas.draw() 
        canvas.get_tk_widget().pack(side=tk.TOP) 
      
        # creates the matplotlib default toolbar 
        toolbar = NavigationToolbar2Tk(canvas, self.root) 
        toolbar.update() 
        canvas.get_tk_widget().pack(side=tk.TOP) 


def process_GITT(GITT_data,settings):
    import numpy as np
    from scipy.stats import linregress
    
    # initial data transformation, time as x-axis, voltage as y-axis
    x = GITT_data['time']
    y = GITT_data['volt']
    
    m_AM = float(settings['m_AM'].entry.get())
    M = float(settings['M_AM'].entry.get())
    theocap = float(settings['theocap'].entry.get())
    c0 = float(settings['c0'].entry.get())
    A = float(settings['A'].entry.get())
    scale = float(settings['scale'].entry.get())
    limiter = float(settings['limiter'].entry.get())
    rho = float(settings['rho'].entry.get())
    
    # assigns capacity if it was provided
    # prefers special capacity over pure capacity, if both are provided

    settings['spec_cap'] = False
    settings['cap'] = False

    if 'spec_cap' in GITT_data:
        print(123123)
        x_spec_cap = GITT_data['spec_cap']
        settings['spec_cap'] = True
    elif 'cap' in GITT_data:
        x_cap = GITT_data['cap']
        settings['cap'] = True  
    
    # calculate specific capacity, current cycle, and ion content at every given time
    charge = True
    current_cycle = 0 # charge-discharge cycle number
    ref_cap = 0
    x_ion = []
    y_cycle = []
    # this logic only works with Biologic output!!!
    # i.e., capacity is reset to 0 for every new charge or discharge cycle
    # capacity only ever rises and never decreases
    if settings['cap']:
        x_spec_cap = []
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
                spec_cap = (ref_cap+cap)/m_AM
            else:
                spec_cap = (ref_cap-cap)/m_AM
            x_spec_cap.append(spec_cap)
            x_ion.append(-spec_cap/theocap + c0)
            
    # this logic works with capacity that rises during charge and decreases during discharge
    # i.e., capacity goes up for charge and down for discharges
    elif settings['spec_cap']:
        charge = True
        for i, spec_cap in enumerate(x_spec_cap):
            x_ion.append(-spec_cap/theocap + c0)
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
    y_deriv_cutoff = limiter*np.mean(list(map(abs, y_deriv)))

    # this part detects when the current is applied and removed
    # makes this less dependent on format of GITT data
    current_on = []
    current_off = []
    on_times = []
    off_times = []
    load = False  
    for i, value in enumerate(y_deriv):
        # detects positive derivative jump
        if y_deriv[i] > abs(scale*y_deriv[i-1]) and y_deriv[i] > y_deriv_cutoff and load == False:
            # switches meaning of jump depending of whether the cell is currently being charged or discharged
            if len(current_on) == 0 or y[current_on[-1]] < y[i]:
                current_on.append(i)
                on_times.append(x[i])
            else:
               current_off.append(i)
               off_times.append(x[i])
            load = True
        # detects negative derivative jump, same logic as for positive derivative jump but flipped
        elif y_deriv[i] < -abs(scale*y_deriv[i-1]) and y_deriv[i] < -y_deriv_cutoff and load == True:
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
    
    V_mol = M/rho    
    # F = 96485.3321 # C/mol
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
        if settings['cap'] or settings['spec_cap']:
            # the +1 is a fix to properly process Arbin data
            D_out['ion'].append(x_ion[current_on[i]+1])
            D_out['spec_cap'].append(x_spec_cap[current_on[i]+1])
            D_out['cycle'].append(y_cycle[current_on[i]+1])   
        D_out['time'].append(x[current_on[i]])
        D_out['volt'].append(y[current_on[i]])
        
        # calculation of the diffusion constants
        dE_s = E4-E1
        dE_t = E3-E2
        D = 4/(np.pi*tau) * (m_AM * V_mol/(M*A))**2 * ((dE_s)/(dE_t))**2
        D_out['diff'].append(D) 
    
    # writes diffusion data derived from GITT measurement
    # write_GITT_data(D_out,settings['cap'] or settings['spec_cap'],ion)

    # for plotting purposes only
    
    if True:
        plot = plot_window(GITT_data,GITT_refined,D_out,settings)

    return D_out

'''
Tkinter GUI stuff
'''

class labeled_entry:
     
     def __init__(self, parent_frame, pos, label, init_value, fs=12):
         
         self._frame = tk.Frame(parent_frame)
         self._frame.grid(row=pos,column=0)
         self._label = ttk.Label(self._frame,width=8)#,font=tk.font.Font(size=fs))
         self._label['text'] = label[0] #string.format(label[0])
         self._label.grid(row=0,column=0)
         
         self._altlabel = ttk.Label(self._frame,width=8)#,font=tk.font.Font(size=fs))
         self._altlabel['text'] = " "+label[1]
         self._altlabel.grid(row=0,column=2)
         
         self.value = tk.StringVar()
         self.entry = ttk.Entry(
             self._frame,
             textvariable=self.value,
             width=8,
             justify='right'
             #font=tk.font.Font(size=fs),
         )
         self.entry.insert(0, init_value)

         self.entry.grid(row=0,column=1) 

'''
creates the main window for loading and saving data
'''
class main_window:
    
    # initializes the base window
    def __init__(self):
        
        self.version = "0.5.0"
        self.icon = "" # obsolete, for compatibility
        
        # styling for font
        self.fs_head = 16
        self.fs_norm = 12
        self.fs_note = 8
        
        self.raw_file = None
        self.raw_filename = 'No raw GITT data loaded.'
        self.GITT_data = 0
        self.D_data = 0
                
        self.settings = {}
        
        self.root = create_window("400x560+120+120","GITT Analysis")
        self.root.columnconfigure(0, weight=1)
        self.frame_buttons()
        self.frame_entry_fields()
        
        self.root.mainloop()

        
    # frame holding the buttons for file management and plotting
    def frame_entry_fields(self):
        
        self.frame_entry_fields = tk.Frame(self.root)
        self.frame_entry_fields.grid(row=1,column=0)
        self.frame_entry_fields.columnconfigure(0, weight=1)
        
        self.entries_title = tk.Label(self.frame_entry_fields,
                                      text = "Properties Active Material (AM)")#,
                                      # font = tk.font.Font(size=self.fs_head))
        self.entries_title.grid(row=0,column=0)
        
        self.settings['m_AM'] = labeled_entry(
            self.frame_entry_fields, 
            pos = 1, 
            label = ['m','g'], 
            init_value = 1,
            fs = self.fs_norm)
        
        self.settings['M_AM'] = labeled_entry(
            self.frame_entry_fields, 
            pos = 2, 
            label = ['M','g/mol'], 
            init_value = 13,
            fs = self.fs_norm)
        
        self.settings['rho'] = labeled_entry(
            self.frame_entry_fields, 
            pos = 3, 
            label = ['ρ','g/cm³'], 
            init_value = 55,
            fs = self.fs_norm)
        
        self.settings['theocap'] = labeled_entry(
            self.frame_entry_fields, 
            pos = 4, 
            label = ['max cap.','mAh/g'], 
            init_value = 100,
            fs = self.fs_norm)
        
        self.settings['ion'] = labeled_entry(
            self.frame_entry_fields, 
            pos = 5, 
            label = ['ion',''], 
            init_value = 'Li',
            fs = self.fs_norm)
        
        self.settings['c0'] = labeled_entry(
            self.frame_entry_fields, 
            pos = 6, 
            label = ['c_0 (ion)',''], 
            init_value = 1.0,
            fs = self.fs_norm)
        
        self.entries_title = tk.Label(self.frame_entry_fields,
                                      text = "Properties Measurement")
                                      # font = tk.font.Font(size=self.fs_head))
        self.entries_title.grid(row=7,column=0)
        
        self.settings['A'] = labeled_entry(
            self.frame_entry_fields, 
            pos = 8, 
            label = ['A_cont','cm²'], 
            init_value = 1.25,
            fs = self.fs_norm)
        
        self.entries_title = tk.Label(self.frame_entry_fields,
                                      text = "Settings Analysis")
                                      # font = tk.font.Font(size=self.fs_head))
        self.entries_title.grid(row=9,column=0)
        
        self.settings['scale'] = labeled_entry(
            self.frame_entry_fields, 
            pos = 10, 
            label = ['scale',''], 
            init_value = 1,
            fs = self.fs_norm)
        
        self.settings['limiter'] = labeled_entry(
            self.frame_entry_fields, 
            pos = 11, 
            label = ['limiter',''], 
            init_value = 0.01,
            fs = self.fs_norm)
        
        
        # print(self.entry_distance.entry.get())
        
    # frame holding the buttons
    def frame_buttons(self):
        
        def top_buttons(self):
            self.frame_top_buttons = tk.Frame(self.root)
            self.frame_top_buttons.grid(row=0,column=0,sticky='NSEW')
            self.frame_top_buttons.columnconfigure(0, weight=1)
            self.frame_top_buttons.columnconfigure(1, weight=1)
            self.frame_top_buttons.columnconfigure(2, weight=1)
            
            sep_top = ttk.Separator(self.frame_top_buttons,orient='horizontal')
            sep_top.grid(row=0,column=0,columnspan=3,sticky='EW')
            
            sep_top2 = ttk.Separator(self.frame_top_buttons,orient='horizontal')
            sep_top2.grid(row=3,column=0,columnspan=3,sticky='EW')
            
            self._button_get_GITT_raw = ttk.Button(self.frame_top_buttons, text = 'Open File', command = lambda : get_GITT_raw())
            self._button_get_GITT_raw.grid(row=1,column=0,sticky='NS')
            
            self._button_process_GITT = ttk.Button(self.frame_top_buttons, text = 'Run Analysis', command = lambda : try_process_GITT())
            self._button_process_GITT.grid(row=1,column=1,sticky='NS')
            
            self._button_save_GITT = ttk.Button(self.frame_top_buttons, text = 'Save File', command = lambda : save_GITT())
            self._button_save_GITT.grid(row=1,column=2,sticky='NS')
            
            self.label_GITT_file = tk.Text(self.frame_top_buttons,
                                          # font = tk.font.Font(size=self.fs_note),
                                          height=2
                                          )
            self.label_GITT_file.insert(tk.END,self.raw_filename)
            self.label_GITT_file.grid(row=2,column=0,columnspan=3,pady=10)
        
        def bottom_buttons(self):
        
            self.frame_bottom_buttons = tk.Frame(self.root)
            self.frame_bottom_buttons.grid(row=2,column=0,sticky='EW')
            self.frame_bottom_buttons.columnconfigure(0, weight=1)
            self.frame_bottom_buttons.columnconfigure(1, weight=1)
            
            sep_bottom = ttk.Separator(self.frame_bottom_buttons,orient='horizontal')
            sep_bottom.grid(row=1,column=0,columnspan=2,sticky='EW')
    
            self._button_help = ttk.Button(self.frame_bottom_buttons, text = 'Help', command = lambda : getting_help())
            self._button_help.grid(row=2,column=0)
            
            self._button_about = ttk.Button(self.frame_bottom_buttons, text = 'About', command = lambda : about())
            self._button_about.grid(row=2,column=1)
            
        def get_GITT_raw():
            filetypes = (
                ('data files', '*.csv;*.txt;*.dat'),
                ('All files', '*.*')
            )
            
            self.raw_file = fd.askopenfilename(
                title='Open GITT raw data file',
                initialdir='./',
                filetypes=filetypes)
            if os.path.isfile(self.raw_file) == True:
                self.GITT_data = fetch_GITT_data(self.raw_file)
                self.raw_filename = 'GITT raw data loaded: '+self.raw_file
                self.frame_top_buttons.destroy()
                top_buttons(self)
            else:
                messagebox.showerror("Error in input file!", "Input file could not be found!")
            return
        
        def try_process_GITT():
            if self.GITT_data == 0:
                messagebox.showerror("No GITT data", "No GITT data loaded!")
            else:
                self.D_data = process_GITT(self.GITT_data,self.settings)
        
        def save_GITT():
            if self.GITT_data == 0:
                messagebox.showerror("No GITT data", "No GITT data loaded!")
            else:
                Files = [('CSV File', '*.csv'),
                    ('All Files', '*.*')]
                savefile = fd.asksaveasfile(filetypes = Files, defaultextension = Files)
                write_GITT_data(savefile,self.D_data,self.settings)
        
        def getting_help():
            print('Help')
            return
        
        def about():
            print('About')
            return
        
        top_buttons(self)
        bottom_buttons(self)

'''
function that ensures that the created windows do not become bigger than the screen
'''
def window_size_limiter(avail_wxh,req_wxh,req_offset_xy):

    actual_wxh = [0,0]
    actual_offsets = [0,0]
    
    # check whether window fits on the current screen with and without offsets
    for i in range(len(avail_wxh)):
        if req_wxh[i] > avail_wxh[i]:
            actual_wxh[i] = avail_wxh[i]
            print("Caution, requested window doesn't fit the screen!")
        elif req_wxh[i] + req_offset_xy[i] > avail_wxh[i]:
            actual_wxh[i] = req_wxh[i]
            actual_offsets[i] = avail_wxh[i] - req_wxh[i]
            print("Caution, requested offset would move window off the screen!")
        else:
            actual_wxh[i] = req_wxh[i]
            actual_offsets[i] = req_offset_xy[i]
    
    return actual_wxh,actual_offsets

'''
function that creates a new window
'''
def create_window(dimensions="500x350+100+100", title = "Tkinter Hello World", icon = ""):
   
    w = int(dimensions.split("x")[0])
    h = dimensions.split("x")[1]
    h = int(h.split("+")[0])
    
    offset_x = int(dimensions.split("+")[1])
    offset_y = int(dimensions.split("+")[2])
    
    # initializes the Tk root window
    window = tk.Tk()
    
    # gets screen properties and centers in upper third
    screen_width = window.winfo_screenwidth()
    screen_height = window.winfo_screenheight()
    
    offset_x = int(screen_width/3 - w / 3)
    offset_y = int(screen_height/3 - h / 3)
    
    # makes sure the window stays within bounds
    actual_wxh, actual_offsets = window_size_limiter([screen_width,screen_height],[w,h], [offset_x,offset_y])
    
    # set a title
    window.title(title)
    
    # specify geometry and max and min measurements
    window.geometry(f"{actual_wxh[0]}x{actual_wxh[1]}+{actual_offsets[0]}+{actual_offsets[1]}")
    window.minsize(10,10)
    window.maxsize(screen_width,screen_height)
    if icon != "":
        window.iconbitmap(icon)
    
    return window





if __name__ == '__main__':
    import tkinter as tk
    from tkinter import ttk
    from tkinter import filedialog as fd
    from tkinter import messagebox
    import os
    import matplotlib as mpl
    
    # mpl.use("TkAgg")
    mpl.use("pgf")
    mpl.use("pdf")
    mpl.use("ps")
    mpl.use("svg")
    
    main = main_window()