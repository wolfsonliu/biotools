# ------------------
# Sequence Tailorbird
# ------------------
# Author: Wolfson Liu
# Date: 2017.12.01
# Version: 0.1.2017.12.01
# Description:
#    Used for multiple sequence processing in GUI.
# ------------------

import tkinter as tk

__version__ = 0.1.2017.12.01

def complement(sequence):
    if not isinstance(sequence, str):
        raise ValueError('Sequence should be string.')
    atcg = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
        'a': 't', 't': 'a', 'g': 'c', 'c': 'g'
    }
    result = ''.join([atcg[x] for x in sequence])
    return result


def reverse_complement(sequence):
    if not isinstance(sequence, str):
        raise ValueError('Sequence should be string.')
    atcg = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
        'a': 't', 't': 'a', 'g': 'c', 'c': 'g'
    }
    result = ''.join([atcg[x] for x in sequence][::-1])
    return result


class MyGUI(tk.Frame):
    def __init__(self, master):
        super().__init__(master)
        self.pack()
        self.in_textbox = self._textbox(
            '序列输入框',
            width=40, textwrap=tk.NONE,
            row=0, column=0
        )
        self.exchange_button = tk.Button(
            self,
            text='输入<=输出',
            command=self.button_exchange
        )
        self.exchange_button.grid(row=0, column=1)
        self.out_textbox = self._textbox(
            '序列输出框',
            width=40, textwrap=tk.NONE,
            row=0, column=2
        )
        entry_frame = tk.Frame(self)
        entry_frame.grid(row=1, column=0)
        self.entry = tk.Entry(
            entry_frame, width=30
        )
        self.entry.grid(row=0, column=0, columnspan=2)
        self.entry_buttons = dict()
        self.entry_buttons['添头'] = tk.Button(
            entry_frame, text='添头',
            command=self.entry_button_addhead
        )
        self.entry_buttons['添头'].grid(row=1, column=0)
        self.entry_buttons['加尾'] = tk.Button(
            entry_frame, text='加尾',
            command=self.entry_button_addtail
        )
        self.entry_buttons['加尾'].grid(row=1, column=1)
        self._buttons(row=1, column=1, columnspan=2)


    def _buttons(self, row, column, columnspan):
        self.buttons = tk.Frame(self)
        self.buttons.grid(
            row=row, column=column, columnspan=columnspan
        )
        self.button_dict = dict()
        self.button_dict['反向'] = tk.Button(
            self.buttons, text='反向', command=self.button_reverse
        )
        self.button_dict['反向'].grid(row=0, column=0)
        self.button_dict['互补'] = tk.Button(
            self.buttons, text='互补', command=self.button_complement
        )
        self.button_dict['互补'].grid(row=0, column=1)
        self.button_dict['大写'] = tk.Button(
            self.buttons, text='大写', command=self.button_upper
        )
        self.button_dict['大写'].grid(row=0, column=2)
        self.button_dict['小写'] = tk.Button(
            self.buttons, text='小写', command=self.button_lower
        )
        self.button_dict['小写'].grid(row=0, column=3)
        self.button_dict['反向互补'] = tk.Button(
            self.buttons,
            text='反向互补',
            command=self.button_reverse_complement
        )
        self.button_dict['反向互补'].grid(row=1, column=0)

    def _textbox(self, label, width, textwrap, row, column):
        # make the text box with x and y scrollbar
        frame = tk.Frame(self)
        frame.grid(row=row, column=column)
        thelabel = tk.Label(frame, text=label)
        thelabel.grid(row=0, column=0)
        thetext = tk.Text(
            frame,
            borderwidth=3,
            width=width,
            relief=tk.SUNKEN,
            wrap=textwrap
        )
        thetext.grid(row=1, column=0)
        thescrollx = tk.Scrollbar(
            frame,
            command=thetext.xview,
            orient=tk.HORIZONTAL
        )
        thescrollx.grid(
            row=2, column=0, sticky=tk.W + tk.E + tk.N + tk.S
        )
        thetext.configure(xscrollcommand=thescrollx.set)
        thescrolly = tk.Scrollbar(
            frame,
            command=thetext.yview,
            orient=tk.VERTICAL
        )
        thescrolly.grid(
            row=1, column=1, sticky=tk.W + tk.E + tk.N + tk.S
        )
        thetext.configure(yscrollcommand=thescrolly.set)
        return thetext

    def button_exchange(self):
        # button exchange function
        output_text = self.out_textbox.get(1.0, tk.END).rstrip()
        self.in_textbox.delete(1.0, tk.END)
        self.in_textbox.insert(1.0, output_text)

    def entry_button_addhead(self):
        entry_text = self.entry.get()
        input_text = self.in_textbox.get(
            1.0, tk.END
        ).rstrip().split('\n')
        output_text = [entry_text + x for x in input_text]
        self.out_textbox.delete(1.0, tk.END)
        self.out_textbox.insert(1.0, '\n'.join(output_text))

    def entry_button_addtail(self):
        entry_text = self.entry.get()
        input_text = self.in_textbox.get(
            1.0, tk.END
        ).rstrip().split('\n')
        output_text = [x + entry_text for x in input_text]
        self.out_textbox.delete(1.0, tk.END)
        self.out_textbox.insert(1.0, '\n'.join(output_text))

    def button_upper(self):
        # button upper function
        output_text = self.in_textbox.get(
            1.0, tk.END
        ).rstrip().upper()
        self.out_textbox.delete(1.0, tk.END)
        self.out_textbox.insert(1.0, output_text)

    def button_lower(self):
        # button lower function
        output_text = self.in_textbox.get(
            1.0, tk.END
        ).rstrip().lower()
        self.out_textbox.delete(1.0, tk.END)
        self.out_textbox.insert(1.0, output_text)

    def button_reverse(self):
        # button reverse function
        input_text = self.in_textbox.get(
            1.0, tk.END
        ).rstrip().split('\n')
        output_text = [x[::-1] for x in input_text]
        self.out_textbox.delete(1.0, tk.END)
        self.out_textbox.insert(1.0, '\n'.join(output_text))

    def button_complement(self):
        # button complement function
        input_text = self.in_textbox.get(
            1.0, tk.END
        ).rstrip().split('\n')
        output_text = [complement(x) for x in input_text]
        self.out_textbox.delete(1.0, tk.END)
        self.out_textbox.insert(1.0, '\n'.join(output_text))

    def button_reverse_complement(self):
        # button reverse complement function
        input_text = self.in_textbox.get(
            1.0, tk.END
        ).rstrip().split('\n')
        output_text = [reverse_complement(x) for x in input_text]
        self.out_textbox.delete(1.0, tk.END)
        self.out_textbox.insert(1.0, '\n'.join(output_text))


if __name__ == '__main__':
    root = tk.Tk()
    root.title('>>> 西昆斯缝叶莺 =*= Sequence Tailorbird <<<')
    app = MyGUI(root)
    root.mainloop()

# ------------------
# EOF
# ------------------
