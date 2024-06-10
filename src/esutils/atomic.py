import pandas  as pd

s = 'Fe4+ 2p'

def print_res(data):
    for i, d in data.iterrows():
        print('{} {} {}'.format(d.E, d.C, d.V))

        # SOC
        print(d.Q.replace('_', ' ').replace('=', ' '))
        if str(d.W) != 'nan':
            print(d.W.replace('_', ' ').replace('=', ' '))
        
        # Fk, Gk
        for j in ['Ie', 'If', 'Ig', 'Ia', 'Ic', 'Ib', 'Id', 'Ih', 'Ii']:
            v = d[j]
            if str(v) != 'nan':
                print(v.replace('=', ' '))
        print('')

def search_data(s):
    names = ['Z', 'C', 'V', 'E', 'Q', 'W', 'Ia', 'Ib', 'Ic', 'Id', 'Ie', 'If', 'Ig', 'Ih', 'Ii']
    data = pd.read_csv('atomic.txt', names=names)
    
    d = 1
    for i in s.split(' '):
        f = 0
        for j in ['C', 'V', 'E']:
            f = f | data[j].str.match(i)
        d = d & f
            
    
    data = data.loc[d, :]
    
    return data
    


data = search_data(s)
print_res(data)