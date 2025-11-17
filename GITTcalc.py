# -*- coding: utf-8 -*-
'''
Created on Sat Feb  8 15:51:35 2025

@author: mhaefner-chem
'''

# INPUT

'''
function reads in the GITT data from a given file
Format: 
    dictionary 'data', keys from list ['time','volt','cap','spec_cap']
    time (continuous measurement time) and volt (measured voltage) are required
    either cap (capacity) or spec_cap (special capacity) are requied
    
    datapoints are stored in lists:
    data[key] = []
'''
def get_GITT_data(file):
    
    data = {}    
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
                    if item.lower() in ['test time (s)','time/s','time (s)','time']:
                        label = 'time'
                    elif item.lower() in ['voltage (v)','ewe/v','volt','voltage']:
                        label = 'volt'
                    elif item.lower() in ['capacity (mah)','capacity/ma.h','capacity/mah','capacity','cap']:
                        label = 'cap'
                    elif item.lower() in ['specific capacity (mah/g)','specificcapacity/ma.h/g','specific_capacity','spec_cap']:
                        label = 'spec_cap'
                        
                    data[label] = []
                    labels.append(label)
            else:
                for item_index, item in enumerate(line.split(splitter)):
                    if labels[item_index] in required:
                        try:
                            float(item)
                        except:
                            messagebox.showerror('Faulty GITT data', 'GITT data contains non-numerical values. Please check the input file.')
                            return 0
                        else:
                            data[labels[item_index]].append(float(item))
                                
                    
    for label in ['time','volt']:
        if not label in data.keys():
            messagebox.showerror('GITT data incomplete', 'The GITT data needs to contain at least a column labeled \'time/s\' and a column labeled \'Ewe/V\'.')
            return 0

    return data

# OUTPUT

'''
This function writes the GITT settings to an .info-file in plain text format associated to the 
filename of the input raw data file. That allows the program to keep the settings between closing
and reopening.
'''
def write_GITT_settings(settings_file, settings):
    
    number_settings = ['refcap','m_AM','M_AM','rho','c0','A','scale','limiter']
    with open(settings_file, mode='w') as f:
        for item in number_settings:  
            value = settings[item].entry_main.get()

            if not item in blacklist_param:
                value += ','
                value += settings[item].entry_error.get()

            f.write('{},{}\n'.format(item,value))
            
'''
This function writes an example file with mock GITT data.
'''
def write_GITT_example():
    import numpy as np
    
    volt_init = 2
    
    t_all = list(np.linspace(0,540,10))
    volt_all = [2]*10
    cap_all = [0]*10
    
    def jump(t_in,volt_in,cap_in,charge):
        t_out = np.linspace(0,8100,136)
        volt_out = []
        cap_out = []
        
        for t in t_out:
            if t < 900:
                volt_out.append((np.sqrt(t)/900 + 0.05)*charge + volt_in)
                cap_out.append(cap_in + t/4500)
            else:
                volt_out.append((np.exp(-(t-300)/900)/20 + 0.03)*charge + volt_in)
                cap_out.append(cap_out[-1])
                
        for i in range(len(t_out)):
            t_out[i] = t_out[i] + t_in
        return t_out, volt_out, cap_out
    
    t_jump, volt_jump, cap_jump = jump(t_all[-1], volt_init, cap_all[-1], 1)
    for i, t in enumerate(t_jump):
        t_all.append(t)
        volt_all.append(volt_jump[i])
        cap_all.append(cap_jump[i])
    
    for i in range(12):
        t_jump, volt_jump, cap_jump = jump(t_all[-1], volt_all[-1], cap_all[-1], 1)
        for i, t in enumerate(t_jump):
            t_all.append(t)
            volt_all.append(volt_jump[i])
            cap_all.append(cap_jump[i])
        
    cap_all[-1] = 0
            
    for i in range(12):
        t_jump, volt_jump, cap_jump = jump(t_all[-1], volt_all[-1], cap_all[-1], -1)
        for i, t in enumerate(t_jump):
            t_all.append(t)
            volt_all.append(volt_jump[i])
            cap_all.append(cap_jump[i])
    
    Files = [('TXT File', '*.txt'),
        ('All Files', '*.*')]
    savefile = fd.asksaveasfile(filetypes = Files, defaultextension = Files)
    
    with savefile as fw:
        fw.write('time/s\tEwe/V\tCapacity/mA.h\n')
        for i, t in enumerate(t_all):
            fw.write('{:.8E}\t{:.8E}\t{:.8E}\n'.format(t,volt_all[i],cap_all[i]))

'''
this function saves the raw and processed GITT data into a OriginLab workbook,
carrying over the name of the original file as label
'''
def write_GITT_2_origin(data_raw,data_out,settings):
    import originpro as op
    book = op.new_book(lname=settings['name'])
    
    # worksheet for raw data
    wks_raw = book.add_sheet(name='raw_data')
    wks_raw.cols = 3
    
    wks_raw.from_list(0,data_raw['time'],lname='Time',units='s')
    wks_raw.from_list(1,data_raw['volt'],lname='Voltage',units='V')
    if settings['cap']:
        wks_raw.from_list(2,data_raw['cap'],lname='Capacity',units='mAh')
    elif settings['spec_cap']:
        wks_raw.from_list(2,data_raw['spec_cap'],lname='Specific Capacity',units='mAh/g')
    
    # worksheet for processed data
    wks_diff = book.add_sheet(name='diffusion_data')
    if settings['cap'] or settings['spec_cap']:
        cols = 7
        half_cycle_label = []
        for value in data_out['cycle']:
            half_cycle_label.append(value + 1)
    else:
        cols = 4
    wks_diff.cols = cols
    
    clm_info = [
        {'data':data_out['time'],'label':'Time','units':'s','comments':''},
        {'data':data_out['volt'],'label':'Voltage','units':'V','comments':''},
        {'data':[x[0] for x in data_out['diff']],'label':'Diffusion Coefficient','units':'cm²/s','comments':''}, 
        {'data':[x[1] for x in data_out['diff']],'label':'ERROR Diffusion Coefficient','units':'cm²/s','comments':''}]
    if settings['cap'] or settings['spec_cap']:
        clm_info.append({'data':[x[0] for x in data_out['spec_cap']],'label':'Specific Capacity','units':'mAh/g','comments':''})
        clm_info.append({'data':[x[0] for x in data_out['ion']],'label':'Content Conducting Ion','units':'','comments':''})
        clm_info.append({'data':half_cycle_label,'label':'Half Cycle','units':'','comments':'odd cycles: charge, even cycles: discharge'})
    
    for idx in range(cols):    
        wks_diff.from_list(idx,
                           clm_info[idx]['data'],
                           lname=clm_info[idx]['label'],
                           units=clm_info[idx]['units'],
                           comments=clm_info[idx]['comments'])
        

    if settings['cap'] or settings['spec_cap']:
        # half cycle worksheets
        current_cycle = -1
        block_idx = []
        for idx, cycle in enumerate(data_out['cycle']):
            if cycle > current_cycle:
                block_idx.append(idx)
                current_cycle = cycle

        wks_cycs = []
        for idx, value in enumerate(block_idx):
            if idx%2 == 0:
                name = 'Charge {}'.format(1+idx//2)
            else:
                name = 'Discharge {}'.format(1+idx//2)
                
            wks_cycs.append(book.add_sheet(name=name))
            wks_cycs[-1].cols = cols-1
            
            if idx == len(block_idx)-1:
                lo_bound = block_idx[idx]
                hi_bound = len(data_out['cycle']) 
            else:
                lo_bound = block_idx[idx]
                hi_bound = block_idx[idx+1]
              
            for iidx in range(cols-1):
                wks_cycs[-1].from_list(iidx,
                                   clm_info[iidx]['data'][lo_bound:hi_bound],
                                   lname=clm_info[iidx]['label'],
                                   units=clm_info[iidx]['units'],
                                   comments=clm_info[iidx]['comments'])           

'''
function outputs diffusion data
Format: 
    dictionary 'data_out', keys from list ['time','volt','spec_cap','cycle','diff']
    'time' is measurement time and 'volt' is measured voltage associated with diffusion rate 'diff'
    
    if capacity data is supplied, output also contains
    'spec_cap' (special capacity)
    'cycle' (current cycle, usually 0 is first charge, 1 is first discharge, 2 is second charge, etc.)
    
    datapoints are stored in lists:
    data_out[key] = []
'''
def write_GITT_data(savefile,data_out,settings,refinement_param):
    
    with savefile as f:
        if settings['cap'] or settings['spec_cap']:            
            f.write('{:16},{:16},{:16},{:16},{:16},{:16},{:16}\n'.format('time/s','volt/V','x_ion','SpecCap/mAh/g','D/cm^2/s','ERR_D/cm^2/s','Cycle'))
            for idx, value in enumerate(data_out['diff']):
                f.write(
                    '{:16.10e},{:16.10e},{:16.10e},{:16.10e},{:16.10e},{:16.10e},{:4}\n'.format(
                        data_out['time'][idx],
                        data_out['volt'][idx],
                        data_out['ion'][idx][0],
                        data_out['spec_cap'][idx][0],
                        data_out['diff'][idx][0],
                        data_out['diff'][idx][1],
                        data_out['cycle'][idx]))
        else:
            f.write('{:16},{:16},{:16},{:16}\n'.format('time/s','volt/V','D/cm^2/s','ERR_D/cm^2/s'))
            for idx, value in enumerate(data_out['diff']):
                f.write('{:16.10e},{:16.10e},{:16.10e},{:16.10e}\n'.format(
                    data_out['time'][idx],
                    data_out['volt'][idx],
                    data_out['diff'][idx][0],
                    data_out['diff'][idx][1]))
    
    # TODO fix or replace
    with open('abc.csv','w') as f:
        print(os.getcwd())
        f.write('E1,D_E1,E2,D_E2,E3,D_E3,E4,D_E4,tau,t_on,t_on,t_off,R^2,t_relax\n')
        for idx, item in enumerate(refinement_param):
            for value in item:
                try: 
                    len(value)
                    f.write('{},{},'.format(value[0],value[1]))
                except:
                    f.write('{},'.format(value))
            f.write('\n')

# METHODS
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
This function processes the raw GITT data
'''
def process_GITT(GITT_data,settings):
    import numpy as np
    from scipy.stats import linregress
    
    # initial data transformation, time as x-axis, voltage as y-axis
    x = GITT_data['time']
    y = GITT_data['volt']
    
    p_list = ['A','m_AM','M_AM','refcap','c0','rho','scale','limiter']
    p_val = {}
    
    for p_key in p_list:
        p_val[p_key] = [0,0]
        try:
            p_val[p_key][0] = float(settings[p_key].entry_main.get())
        except:
            messagebox.error('Faulty Settings','The setting {} contains non-numerical data. Please check setting and correct.'.format(p_key))
            return 0,0
    
    # assigns capacity if it was provided
    # prefers special capacity over pure capacity, if both are provided
    settings['spec_cap'] = False
    settings['cap'] = False

    if 'spec_cap' in GITT_data:
        raw_spec_cap = GITT_data['spec_cap']
        settings['spec_cap'] = True
    elif 'cap' in GITT_data:
        x_cap = GITT_data['cap']
        settings['cap'] = True  
    
    p_val['m_AM'][0] = p_val['m_AM'][0]/ 1000
    p_val['m_AM'][1] = p_val['m_AM'][1]/ 1000
    
    def calculate_V_mol(M_AM,rho):
        
        V_mol = M_AM[0] / rho[0]
        
        d_M_AM = M_AM[1]/(rho[0])
        d_rho = rho[1]*M_AM[0]/(rho[0]**2)        
        V_mol_err = np.sqrt(d_M_AM**2+d_rho**2)
        
        return V_mol, V_mol_err
    
    p_val['V_mol'] = calculate_V_mol(p_val['M_AM'],p_val['rho'])
    
    
    
    # calculate specific capacity, current cycle, and ion content at every given time
    charge = True
    current_cycle = 0 # charge-discharge cycle number
    ref_cap = 0
    max_cap = 0
    x_ion = []
    x_spec_cap = []
    y_cycle = []    
    
    # this logic works with regular output from a BioLogic cycler
    # i.e., capacity is reset to 0 for every new charge or discharge cycle
    # capacity only ever rises and never decreases
    def calculate_x_ion(spec_cap,ref_cap,c0):
        x_ion_value = -spec_cap[0]/ref_cap[0] + c0[0]
                       
        d_spec_cap = spec_cap[1]/ref_cap[0]
        d_ref_cap = ref_cap[1]*spec_cap[0]/ref_cap[0]**2
        d_c0 = c0[1]
        x_ion_err = np.sqrt(d_spec_cap**2+d_ref_cap**2+d_c0**2)
        
        return [x_ion_value, x_ion_err]
    
    if settings['cap']:
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
            spec_cap = [0,0]
            if charge == True:
                spec_cap[0] = (ref_cap+cap)/p_val['m_AM'][0]
            else:
                spec_cap[0] = (ref_cap-cap)/p_val['m_AM'][0]
            spec_cap[1] = p_val['m_AM'][1]/p_val['m_AM'][0]**2
            x_spec_cap.append(spec_cap)
            
            x_ion.append(calculate_x_ion(spec_cap,p_val['refcap'],p_val['c0']))
            
            
    elif settings['spec_cap']:
        charge = True

        for i in range(len(raw_spec_cap)):
            y_cycle.append(current_cycle)
            if i > 1:
                if raw_spec_cap[i] > max_cap:
                    max_cap = raw_spec_cap[i]
                    
                if raw_spec_cap[i] > 1E-4 and raw_spec_cap[i] < 1E-4*max_cap:
                    current_cycle += 1
                    if charge:
                        charge = False
                        ref_cap = ref_cap+max_cap  
                    else:
                        charge = True
                        ref_cap = ref_cap-max_cap
                    max_cap = 0
                
                spec_cap = [0,0]
                if charge == True:
                    spec_cap[0] = (ref_cap+raw_spec_cap[i])
                else:
                    spec_cap[0] = (ref_cap-raw_spec_cap[i])
                    
                x_spec_cap.append(spec_cap) # cannot compute error without errors from measurement yet
                x_ion.append(calculate_x_ion(spec_cap,p_val['refcap'],p_val['c0']))

    # get numerical derivative of voltage
    # cutoff determines minimum jump in derivative required for it to be counted
    y_deriv = get_numerical_derivative(x, y)
    y_deriv_cutoff = p_val['limiter'][0]*np.mean(list(map(abs, y_deriv)))

    # this part detects when the current is applied and removed
    # makes this less dependent on format of GITT data
    current_on = []
    current_off = []
    on_times = []
    off_times = []
    load = False
    bad_fit = 0
    bad_expol = 0    
    n_datapoints = len(y_deriv)
    
    for i, value in enumerate(y_deriv):
        
        # section to make sure that the new y-value in the comparison is 60 s after the voltage jump after current application
        # ensures that charge/discharge is recognized more accurately
        delta_i = 1
        delta_t = 0
        while delta_t < 60:
            if delta_i + i < n_datapoints:
                delta_t = x[delta_i+i] - x[i]
                delta_i += 1
            else:
                break
            
        # detects positive derivative jump
        if y_deriv[i] > abs(p_val['scale'][0]*y_deriv[i-1]) and y_deriv[i] > y_deriv_cutoff and load == False:
            # switches meaning of jump depending of whether the cell is currently being charged or discharged
            if len(current_on) == 0 or y[current_on[-1]] < y[i+delta_i]: #y_deriv[i] > 0: #
                current_on.append(i)
                on_times.append(x[i])
            else:
               current_off.append(i)
               off_times.append(x[i])
            load = True
        # detects negative derivative jump, same logic as for positive derivative jump but flipped
        elif y_deriv[i] < -abs(p_val['scale'][0]*y_deriv[i-1]) and y_deriv[i] < -y_deriv_cutoff and load == True:
            if len(current_on) == 0 or y[current_on[-1]] < y[i+delta_i]: #y_deriv[i] < 0: #
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
    GITT_refined = []
    
    off = 0
    for i, on in enumerate(current_on):   
        # discard the first cycle
        if i == 0:
            continue
        
        # determines the next point after on at which current is turned off
        relax = x[on]-x[off]
        for off in current_off:
            if off > on:
                break
        
        # determination of E1, E3, and tau
        E1 = [0,0] # currently no error evaluated
        E1[0] = y[on]
        tau = x[off]-x[on]
        
        # E4 requires logic to properly treat final titration
        E4 = [0,0] # currently no error evaluated
        if i < len(current_on)-1:
            E4[0] = y[current_on[i+1]]
        else:
            E4[0] = y[-1]
            
        # E2 requires linear regression for sqrt-behavior while current is applied
        interval_x = []
        interval_y = []
        test_x = []
        for j in range(on,off):
            if x[j]-x[on] > tau/5:
                interval_x.append(np.sqrt(x[j]-x[on]))
                test_x.append(x[j]-x[on])
                interval_y.append(y[j])
        try:
            regress_param = linregress(interval_x,interval_y)
        except:
            continue
        if any(np.isnan(regress_param)):
            continue
        
        m = [regress_param.slope,regress_param.stderr]
        b = [regress_param.intercept,regress_param.intercept_stderr]
        
        E2 = [0,0]
        E2[0] = b[0]
        E2[1] = b[1]
        
        if regress_param.rvalue**2 < 0.99:
            bad_fit += 1
        
        E3 = [0,0]
        E3[0] = m[0] * np.sqrt(tau) + b[0]
        E3[1] = np.sqrt((m[1]*np.sqrt(tau))**2 + b[1]**2)
        
        # check if E2 is lower than E1 during charge or higher during discharge, set E2 to E1
        if E2[0] < E3[0] and E2[0] < E1[0]:
            E2[0] = E1[0]
            E2 = E1
            bad_expol += 1
        elif E2[0] > E3[0] and E2[0] > E1[0]:
            E2 = E1
            bad_expol += 1

        
        GITT_refined.append((E1,E2,E3,E4,tau,x[on],x[on],x[off],regress_param.rvalue**2,relax)) # required for plotting
        
        # collect data for output
        if settings['cap'] or settings['spec_cap']:
            # the +1 is a fix to properly process Arbin data
            D_out['ion'].append(x_ion[current_on[i]+1])
            D_out['spec_cap'].append(x_spec_cap[current_on[i]+1])
            D_out['cycle'].append(y_cycle[current_on[i]+1])   
        D_out['time'].append(x[current_on[i]])
        D_out['volt'].append(y[current_on[i]])
        
        # calculation of the diffusion constants
        dE_s = E4[0]-E1[0]
        dE_t = E3[0]-E2[0]
        
        D = [0,0]
        D[0] = 4/(np.pi*tau) * (p_val['m_AM'][0] * p_val['V_mol'][0]/(p_val['M_AM'][0]*p_val['A'][0]))**2 * ((dE_s)/(dE_t))**2
        
        # errors: m_AM, V_mol, M_AM, A, E2, E3
        
        pf = 8/(np.pi*tau)
        d_m_AM = 2*p_val['m_AM'][0] * pf*(p_val['V_mol'][0]/(p_val['A'][0]*p_val['M_AM'][0]))**2 * p_val['A'][0]**-3*((dE_s)/(dE_t))**2
        d_m_AM *= p_val['m_AM'][1]
        
        d_V_mol = 2*p_val['V_mol'][0] * pf*(p_val['m_AM'][0]/(p_val['A'][0]*p_val['M_AM'][0]))**2 * p_val['A'][0]**-3*((dE_s)/(dE_t))**2
        d_V_mol *= p_val['V_mol'][1]
        
        d_M_AM = pf*(p_val['m_AM'][0]*p_val['V_mol'][0]/p_val['A'][0])**2 * p_val['M_AM'][0]**-3*((dE_s)/(dE_t))**2
        d_M_AM *= p_val['M_AM'][1]
        
        
        d_A = pf*(p_val['m_AM'][0]*p_val['V_mol'][0]/p_val['M_AM'][0])**2 * p_val['A'][0]**-3*((dE_s)/(dE_t))**2
        d_A *= p_val['A'][1]
        
        tmp = pf*(p_val['m_AM'][0]*p_val['V_mol'][0]/(p_val['M_AM'][0]*p_val['A'][0]))**2*((dE_s)/(dE_t))**2/dE_t
        d_E2 = tmp*E2[1]
        d_E3 = tmp*E3[1]

        
        D[1] = np.sqrt(d_m_AM**2+d_V_mol**2+d_M_AM**2+d_A**2+d_E2**2+d_E3**2)
        
        
        D_out['diff'].append(D) 
    
    # check whether there is issues with the titration lengths
    def evaluate_tau(taus,relaxes):
        
        buckets = [0,0,0]
        tot_length = [0,0,0]
        median_tau = np.median(taus)
        median_rlx = np.median(relaxes)
        
        for tau in taus:
            if tau > 0.95*median_tau and tau < 1.05*median_tau:
                buckets[1] += 1
                tot_length[1] += tau
            elif tau > 1.05*median_tau:   
                buckets[2] += 1
                tot_length[2] += tau
            elif tau < 0.95*median_tau:
                buckets[0] += 1
                tot_length[0] += tau
        
        if tot_length[2] > tot_length[1] or buckets[2] > 0.01*buckets[1]:
            messagebox.showwarning('Check Results','A significant number of abnormally long titration cycles was obtained, indicating that the program failed to correctly identify all titration cycles. Please reduce the settings \'scale\' and \'limiter\'.')
        
        if median_tau > median_rlx/4:
            messagebox.showwarning('Titration timings','''Titration duration
Current applied for {:.2f} s
Relaxing for {:.2f} s
Relaxation time is relatively short compared to charging time. Make sure that the relaxation is sufficiently long!'''.format(median_tau,median_rlx))
        elif settings['timing'].get():
            messagebox.showinfo('Titration timings','''Titration duration
Current applied for {:.2f} s
Relaxing for {:.2f} s'''.format(median_tau,median_rlx))

    evaluate_tau(list(zip(*GITT_refined))[4],list(zip(*GITT_refined))[9])
    
    if bad_fit > 0:
        messagebox.showinfo('Check Results','The regression for determining the onset energy yielded a bad fit {} times. Please check the results for errors and outliers.'.format(bad_fit))

    return D_out, GITT_refined

# GUI

'''
Function for plotting results from analyzed GITT data
messy code, refactoring not yet planned
'''
class plot_window:
    # initializes the window and default plotting data
    def __init__(self,GITT_data,GITT_refined,D_out,settings):
        
        self.root = create_window('1000x700+120+120', 'V-t and D-t plot')
    
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
        from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,  
        NavigationToolbar2Tk)
        import matplotlib as mpl
        import numpy as np
        from math import log10, floor
        
        plt.close('all')
        mpl.use('pgf')
        mpl.use('pdf')
        mpl.use('ps')
        mpl.use('svg')
        
        fs = 24
        
        fig = plt.figure(figsize=(2.5,2.5), dpi=300)
        ax = fig.add_subplot()
        
        colors = ['green','blue','red']
        alphas = [1.0,1.0,1.0,1.0]
        labels = ['E1','E2','E3']
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        for i in range(3):
            tmp = [[],[],[]]
            for result in self.refined:
                if i < 3:
                    tmp[0].append(result[i+5])
                    tmp[1].append(result[i][0])
                    tmp[2].append(result[i][1])
                else:
                    tmp[0].append(result[i+4])
                    tmp[1].append(result[-1])
                    tmp[2].append(0)
                    
            ax.scatter(tmp[0],tmp[1],marker='x',color=colors[i],s=100,zorder=50,label=labels[i],alpha=alphas[i])
            ax.errorbar(tmp[0],tmp[1],xerr=0,yerr=tmp[2],fmt='none',color=colors[i])
        
        for result in self.refined:
            if result[-2] < 0.99:
                ax.text(result[6],result[1][0],'  R² {:6.4}'.format(result[-2]),
                        ha='left', va='center',
                        fontsize=fs,bbox=props)
        
          
        ax.plot(self.data['time'],self.data['volt'],linestyle='-',label='E',marker='x',markersize=1)
        
        times = self.D['time']
        x_max_tick = round(max(times), -int(floor(log10(max(times)))))
        plt.xticks(np.arange(0, x_max_tick, x_max_tick/5),fontsize=fs)
        plt.yticks(fontsize=fs)
        
        ax2 = ax.twinx()
        ax2.scatter(self.D['time'],[x[0] for x in self.D['diff']],color='black',marker='+',label='D')
        ax2.errorbar(self.D['time'],[x[0] for x in self.D['diff']],xerr=0,yerr=[x[1] for x in self.D['diff']],fmt='none',color='black')
        
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax2.legend(h1+h2, l1+l2,loc=1,fontsize=fs)
        
        ax2.set_yscale('log')
        
        filtered_D = np.ma.masked_invalid(self.D['diff'])
        ylim2_max = max([x[0] for x in filtered_D])*5
        ylim2_min = min([x[0] for x in filtered_D])/5
        ax2.set_ylim(ylim2_min,ylim2_max)
        
        plt.yticks(fontsize=fs)
        plt.tight_layout()
        
        
        ax.xaxis.set_major_formatter(mpl.ticker.StrMethodFormatter("{x:.0f}"))

        ax.grid(zorder=-50,linestyle='--',alpha=0.66)
        ax.set_xlabel('Time t (s)',fontsize=fs) 
        ax.set_ylabel('Voltage E (V)',fontsize=fs)  
        ax2.set_ylabel('Diffusion rate D (cm²/s)',fontsize=fs) 
        ax.set_title('GITT Plot\nDiffusion and Voltage against Time',fontsize=fs+4)
    
        # creates and places Tkinter canvas for the matplotlib figure
        canvas = FigureCanvasTkAgg(fig, master = self.root)   
        canvas.draw() 
        canvas.get_tk_widget().pack(side=tk.TOP,fill='both',expand=False)
        
        # creates the matplotlib default toolbar 
        toolbar = NavigationToolbar2Tk(canvas, self.root) 
        toolbar.update()
        canvas.get_tk_widget().pack(side=tk.TOP,fill='both',expand=True) 
      
'''
Function to streamline the creation of the labeled entries
'''
class labeled_entry:
     
     def __init__(self, parent_frame, pos, label, init_value, b_error=False):
         
         self._frame = tk.Frame(parent_frame)
         self._frame.grid(row=pos,column=0)
         
         self._label = ttk.Label(self._frame,width=8)
         self._label['text'] = label[0]
         self._label.grid(row=0,column=0)
         
         self.value = [tk.StringVar(),tk.StringVar()]
         
         self.entry_main = ttk.Entry(
             self._frame,
             textvariable=self.value[0],
             width=12,
             justify='right'
         )
         self.entry_main.insert(0, init_value)
         self.entry_main.grid(row=0,column=1) 
         
         self._altlabel = ttk.Label(self._frame,width=8)
         self._altlabel['text'] = ' '+label[1]
         
         if b_error:
             self._label_error = ttk.Label(self._frame)
             self._label_error['text'] = '+/-'
             self._label_error.grid(row=0,column=2)
             
             self.entry_error = ttk.Entry(
                 self._frame,
                 textvariable=self.value[1],
                 width=6,
                 justify='right'
             )
             self.entry_error.insert(0, 0)
             self.entry_error.grid(row=0,column=3) 
             
             self._altlabel.grid(row=0,column=4)
         else:
            self._altlabel.grid(row=0,column=2)

         

'''
Creates the main window for loading and saving data
'''
class main_window:
    
    # initializes the base window
    def __init__(self):
        
        self.version = '0.9.8'
        self.icon = ''
        
        self.raw_file = None
        self.raw_filename = ''
        self.GITT_data = 0
        self.D_data = 0
                
        self.settings = {}
        
        self.root = create_window('450x550+120+120','GITTcalc')
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(1, weight=1)
        self.frame_buttons()
        self.frame_entry_fields()
        
        self.root.mainloop()
  
    '''
    This frame contains all entry fields for the settings
    '''
    def frame_entry_fields(self):
        
        self.frame_entry_fields = tk.Frame(self.root)
        self.frame_entry_fields.grid(row=1,column=0,sticky='N')
        self.frame_entry_fields.columnconfigure(0, weight=1)
        
        self.entries_title = tk.Label(self.frame_entry_fields,
                                      text = 'Properties Active Material (AM)'
                                      )
        self.entries_title.grid(row=1,column=0)
        
        
        
        self.settings['m_AM'] = labeled_entry(
            self.frame_entry_fields, 
            pos = 1, 
            label = ['m','mg'], 
            init_value = 5,
            b_error = True)
        
        self.settings['M_AM'] = labeled_entry(
            self.frame_entry_fields, 
            pos = 2, 
            label = ['M','g/mol'], 
            init_value = 100,
            b_error = False)
        
        self.settings['rho'] = labeled_entry(
            self.frame_entry_fields, 
            pos = 3, 
            label = ['ρ','g/cm³'], 
            init_value = 4,
            b_error = True)
        
        self.settings['refcap'] = labeled_entry(
            self.frame_entry_fields, 
            pos = 4, 
            label = ['ref cap.','mAh/g'], 
            init_value = 150,
            b_error = False)
        
        self.settings['c0'] = labeled_entry(
            self.frame_entry_fields, 
            pos = 5, 
            label = ['c_0 (ion)',''], 
            init_value = 1.0,
            b_error = False)
        
        self.entries_title = tk.Label(self.frame_entry_fields,
                                      text = 'Properties Measurement')
        self.entries_title.grid(row=6,column=0)
        
        self.settings['A'] = labeled_entry(
            self.frame_entry_fields, 
            pos = 7, 
            label = ['A_cont','cm²'], 
            init_value = 1.25,
            b_error = True)
        
        self.entries_title = tk.Label(self.frame_entry_fields,
                                      text = 'Settings Analysis')
        self.entries_title.grid(row=8,column=0)
        
        self.settings['scale'] = labeled_entry(
            self.frame_entry_fields, 
            pos = 9, 
            label = ['scale',''], 
            init_value = 1)
        
        self.settings['limiter'] = labeled_entry(
            self.frame_entry_fields, 
            pos = 10, 
            label = ['limiter',''], 
            init_value = 0.01)
        
        self.frame_checkbts = tk.Frame(self.frame_entry_fields)
        self.frame_checkbts.grid(row=11,column=0,sticky='N')
        
        self.settings['plot'] = tk.BooleanVar()
        _checkbt_plot = tk.Checkbutton(self.frame_checkbts, text = "Plot results", 
                variable = self.settings['plot'], 
                onvalue = True, 
                offvalue = False)
        _checkbt_plot.select()
        _checkbt_plot.grid(row=0,column=0,sticky='W')
        
        self.settings['timing'] = tk.BooleanVar()
        _checkbt_plot = tk.Checkbutton(self.frame_checkbts, text = "Print timings", 
                variable = self.settings['timing'], 
                onvalue = True, 
                offvalue = False)
        _checkbt_plot.grid(row=0,column=1,sticky='W')
        
    '''
    This frame handles all relevant buttons.
    '''
    def frame_buttons(self):

        '''
        This is the top row of buttons.
        '''
        def top_buttons(self):
            columns = 3
            
            self.frame_top_buttons = tk.Frame(self.root)
            self.frame_top_buttons.grid(row=0,column=0,sticky='NSEW')
            for i in range(columns):
                self.frame_top_buttons.columnconfigure(i, weight=1)
            
            sep_top = ttk.Separator(self.frame_top_buttons,orient='horizontal')
            sep_top.grid(row=0,column=0,columnspan=columns,sticky='EW')
            
            _button_get_GITT_raw = ttk.Button(self.frame_top_buttons,
                                              text = 'Open File',
                                              command = lambda : get_GITT_raw())
            _button_get_GITT_raw.grid(row=1,column=0,sticky='NS')
            
            _button_process_GITT = ttk.Button(self.frame_top_buttons,
                                              text = 'Run Analysis',
                                              command = lambda : try_process_GITT())
            _button_process_GITT.grid(row=1,column=1,sticky='NS')
            
            _button_save_GITT = ttk.Button(self.frame_top_buttons,
                                           text = 'Save File',
                                           command = lambda : save_GITT())
            _button_save_GITT.grid(row=1,column=2,sticky='NS')
            
            try:
                import originpro as op
                op.org_ver()
                # _button_save_GITT['state'] = 'disabled'
            except:
                pass
            
            sep_top2 = ttk.Separator(self.frame_top_buttons,orient='horizontal')
            sep_top2.grid(row=3,column=0,columnspan=columns,sticky='EW')
            
            
            label_GITT_file = tk.Text(self.frame_top_buttons,
                                          height=4)
            if self.GITT_data == 0:
                label_GITT_file.insert(tk.END,'No raw GITT data loaded.')
            else:
                label_GITT_file.insert(tk.END,self.raw_filename)
            if self.D_data != 0:
                label_GITT_file.insert(tk.END,'\nData processed!')
            
            label_GITT_file.grid(row=2,column=0,columnspan=3,pady=10)
        
        '''
        This is the bottom row of buttons.
        '''
        def bottom_buttons(self):
        
            columns = 2
            self.frame_bottom_buttons = tk.Frame(self.root)
            self.frame_bottom_buttons.grid(row=2,column=0,sticky='EW')
            for i in range(columns):
                self.frame_bottom_buttons.columnconfigure(i, weight=1)
            
            sep_bottom = ttk.Separator(self.frame_bottom_buttons,orient='horizontal')
            sep_bottom.grid(row=1,column=0,columnspan=columns,sticky='EW')
    
            _button_help = ttk.Button(self.frame_bottom_buttons,
                                      text = 'Help',
                                      command = lambda : getting_help())
            _button_help.grid(row=2,column=0,sticky='S')
            
            _button_about = ttk.Button(self.frame_bottom_buttons,
                                       text = 'About',
                                       command = lambda : about())
            _button_about.grid(row=2,column=1,sticky='S') 
        
        '''
        This function retrieves the settings from an existing file.
        '''
        def fetch_GITT_settings():
            
            number_settings = ['refcap','m_AM','M_AM','rho','c0','A','scale','limiter']
            settings_file = self.raw_file.split('.')[0]+'.info'
            if os.path.isfile(settings_file):
                with open(settings_file, mode='r') as f:
                    for line in f.readlines():
                        for item in number_settings:
                            if line.split(',')[0] == item:
                                self.settings[item].entry_main.delete(0,tk.END)
                                self.settings[item].entry_main.insert(0, float(line.split(',')[1]))
                                if not item in blacklist_param:
                                    self.settings[item].entry_error.delete(0,tk.END)
                                    self.settings[item].entry_error.insert(0, float(line.split(',')[2]))
    
        '''
        This function imports raw GITT data.
        '''
        def get_GITT_raw():
            filetypes = (
                ('data files', '*.csv;*.txt;*.dat'),
                ('All files', '*.*'))
            
            self.raw_file = fd.askopenfilename(
                title='Open GITT raw data file',
                initialdir='./',
                filetypes=filetypes)
            
            if os.path.isfile(self.raw_file) == True:
                self.GITT_data = get_GITT_data(self.raw_file)
                if self.GITT_data != 0:
                    self.raw_filename = 'GITT raw data loaded: '+self.raw_file
                    self.frame_top_buttons.destroy()
                    top_buttons(self)
                    fetch_GITT_settings()
                    self.settings['name'] = os.path.basename(self.raw_file.split('.')[0])                
            elif not self.raw_file == '':
                messagebox.showerror('No input file!', 'Input file could not be found!')
        
        '''
        This function processes raw GITT data. Depending on whether launched in an OriginLab
        environment or not, it either writes the results to an OriginLab workbook or
        plots the D-t and V-t diagram for checking the sensibility of the results
        '''
        def try_process_GITT():
            if self.GITT_data == 0:
                messagebox.showerror('No GITT data', 'No GITT data loaded!')
            else:                
                self.D_data, self.GITT_extra = process_GITT(self.GITT_data,self.settings)
                self.frame_top_buttons.destroy()
                top_buttons(self)
                file = self.raw_file.split('.')[0]+'.info'
                write_GITT_settings(file, self.settings)
                try:
                    import originpro as op
                    op.org_ver()
                    has_originpro = True
                except:
                    has_originpro = False
                
                if has_originpro:
                    write_GITT_2_origin(self.GITT_data,self.D_data,self.settings)
                
                if self.settings['plot'].get():
                    plot_window(self.GITT_data,self.GITT_extra,self.D_data,self.settings)
                
        '''
        This function handles GUI side of saving the processed GITT data.
        '''
        def save_GITT():
            if self.GITT_data == 0:
                messagebox.showerror('No GITT data', 'No GITT data loaded!')
                return
            elif self.D_data == 0:
                self.D_data, self.GITT_extra = process_GITT(self.GITT_data,self.settings)
                self.frame_top_buttons.destroy()
                top_buttons(self)
            
            Files = [('CSV File', '*.csv'),
                ('All Files', '*.*')]
            savefile = fd.asksaveasfile(filetypes = Files, defaultextension = Files)
            write_GITT_data(savefile,self.D_data,self.settings,self.GITT_extra)
        
        '''
        This function handles the window containing an overview about the formatting of the raw GITT data input file and the meaning of the different required settings for processing.
        '''
        def getting_help():
            help_frame = create_window('850x550+120+120', 'Quick Tips',self.icon)
            help_frame.config(bg='#FFFFFF')
            message ='''Format raw data:
Natively works with standard output format from BioLogic cyclers. An example raw data file in that format can be generated with the button below. As alternative to tab stops, commas (regular CSV file formatting) or spaces can be used in the input file, e.g.,
0.0E+00 1.0E+00 1.5E+00
0.0E+00,1.0E+00,1.5E+00

Format capacity data:
If capacity should be considered in the analysis, it needs to be formatted such that every half cycle (charge and discharge) are treated separately, each beginning at capacity 0. Examplary behavior for a full cycle (charge and discharge) is simulated in an idealized fashion in the example raw data file.
Regular capacity is indicated by 'Capacity/mA.h' in the input file header,
specific capacity is indicated by 'SpecialCapacity/mA.h/g' in the input file header, e.g.,
Time/s,Ewe/V,Capacity/mA.h,SpecificCapacity/mA.h/g
0.0E+00,1.0E+00,1.5E+00,1.5E+01

Capacity data is not required for analysis. Only time and voltage are required.

Settings for analysis:
m
    mass of the active material in g/cm²
M
    molar mass of the active material in g/mol
ρ
    densitiy of the active material in g/cm³
A_cont
    contact area for measurement in cm²
    
With capacity data:
ref cap.
    specific capacity at ion content 1, e.g., Li1NiO2 or Na1CoO2
c_0
    actual starting ion content

Settings for pattern recognition:
scale (regular range: 1-3)
    decrease if not all titrations are detected
limiter (regular range: 0.02-0.05)
    decrease after scale until all titrations are detected
    adjust until the smoothest curve for diffusion coefficients is obtained
    
Processed data is saved in CSV format.
If program is run in OriginLab either via
    run -pyf GITTcalc.py
or via the OPX file, the results are automatically filled into a new workbook upon analysis.
If program is run as standalone, results are automatically plotted upon analysis.
'''

            text_box = tk.Text(help_frame, wrap = 'word')
            text_box.pack(expand=False,fill=tk.X)
            text_box.insert('end', message)
            text_box.config(state='disabled')
            
            button_example = ttk.Button(help_frame,
                                        text = 'Make Example Input',
                                        command = lambda : write_GITT_example())
            button_example.pack(side=tk.BOTTOM,expand=False,fill=tk.X) 
        
        '''
        This function handles the window containing a short description, version, and license of the program.
        '''
        def about():
            about_frame = create_window('850x600+120+120', 'About BondFinder',self.icon)
            about_frame.config(bg='#FFFFFF')
            message ='''GITTcalc, version {}
            
GITTcalc analyzes raw GITT data to extract the diffusion coefficients of the conducting ion based on the mini-review 'Principle and Applications of Galvanostatic Intermittent Titration Technique for Lithium-ion Batteries' by Jaeyoung Kim, Sangbin Park, Sunhyun Hwang, and Won-Sub Yoon. (DOI: https://doi.org/10.33961/jecst.2021.00836) and equation 16, in particular.

MIT License
Copyright (c) 2025 mhaefner-chem
Contact: michael.haefner@uni-bayreuth.de

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the 'Software'), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''.format(self.version)

            text_box = tk.Text(about_frame, wrap = 'word', height=50)
            text_box.pack(expand=True,fill=tk.X)
            text_box.insert('end', message)
            text_box.config(state='disabled')
    
        top_buttons(self)
        bottom_buttons(self)

'''
function that creates a new window
'''
def create_window(dimensions='500x350+100+100', title = 'You should not be reading this', icon = ''):
   
    w = int(dimensions.split('x')[0])
    h = dimensions.split('x')[1]
    h = int(h.split('+')[0])
    
    offset_x = int(dimensions.split('+')[1])
    offset_y = int(dimensions.split('+')[2])
    
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
    window.geometry(f'{actual_wxh[0]}x{actual_wxh[1]}+{actual_offsets[0]}+{actual_offsets[1]}')
    window.minsize(10,10)
    window.maxsize(screen_width,screen_height)
    if icon != '':
        window.iconbitmap(icon)
    
    return window

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
        elif req_wxh[i] + req_offset_xy[i] > avail_wxh[i]:
            actual_wxh[i] = req_wxh[i]
            actual_offsets[i] = avail_wxh[i] - req_wxh[i]
        else:
            actual_wxh[i] = req_wxh[i]
            actual_offsets[i] = req_offset_xy[i]
    
    return actual_wxh,actual_offsets

if __name__ == '__main__':
    import tkinter as tk
    from tkinter import ttk
    from tkinter import filedialog as fd
    from tkinter import messagebox
    import os
    
    global blacklist_param
    blacklist_param = ['scale','limiter','refcap','c0','M_AM']
    
    main = main_window()