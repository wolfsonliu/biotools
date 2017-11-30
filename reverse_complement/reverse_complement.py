import tkinter as tk

def reverse_complement(sequence):
    if not isinstance(sequence, str):
        raise ValueError('sequence should be string.')
    atcg = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    result = ''.join([atcg[x] for x in sequence.upper()][::-1])
    if sequence.isupper():
        return result
    else:
        return result.lower()

class MyGUI(tk.Frame):
    def __init__(self, master):
        super().__init__(master)
        self.pack()
        self.inputtext = tk.Text(
            master,
            borderwidth=3,
            relief=tk.SUNKEN
        )
        self.inputtext.pack(
            side=tk.LEFT, expand=True, fill=tk.BOTH,
            ipadx=2, ipady=2, padx=2, pady=2
        )
        scrolli = tk.Scrollbar(master, command=self.inputtext.yview)
        scrolli.pack(side=tk.LEFT, expand=True)
        self.inputtext['yscrollcommand'] = scrolli.set
        self.outputtext = tk.Text(
            master,
            borderwidth=3,
            relief=tk.SUNKEN
        )
        self.outputtext.pack(
            side=tk.RIGHT, expand=True, fill=tk.BOTH,
            ipadx=2, ipady=2, padx=2, pady=2
        )
        scrollo = tk.Scrollbar(master, command=self.outputtext.yview)
        scrollo.pack(side=tk.RIGHT, expand=True)
        self.outputtext['yscrollcommand'] = scrollo.set
        self.calculate = tk.Button(
            master,
            text='Reverse Complement',
            command=self.rc
        )
        self.calculate.pack(
            side=tk.BOTTOM, padx=2, pady=2, ipadx=2, ipady=2
        )
    def rc(self):
        self.input = self.inputtext.get(1.0, tk.END).split('\n')
        self.result = [
            reverse_complement(x) for x in self.input
        ]
        self.outputtext.insert(tk.END, '\n'.join(self.result))



if __name__ == '__main__':
    root = tk.Tk()
    app = MyGUI(root)
    root.mainloop()
