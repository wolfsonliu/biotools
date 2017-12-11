# ------------------
# Woodpecker
# ------------------
# Author: Wolfson Liu
# Date: 2017.12.11
# Version: 0.1.2017.12.11
# Description:
#    Used for dowloading multiple papers in GUI.
# ------------------
__version__ = '0.1.2017.12.11'

import os
import bs4
import time
import requests
import tkinter as tk

from requests import ConnectionError

class SciHubError(ValueError):
    pass


class ArticleNotFound(ValueError):
    pass

# ------------------

scihubs = [
    'http://www.sci-hub.tw/',
    'http://www.sci-hub.cc/',
    'http://www.sci-hub.bz/',
    'http://www.sci-hub.hk/',
]


def getscihub(scihubs):
    for sh in scihubs:
        try:
            shpage = requests.get(sh)
            return sh
        except ConnectionError:
            continue
    return sh

goodscihub = getscihub(scihubs)

def scihub_pdfurl(doi, scihub):
    # get scihub pdf url
    shpage = requests.get(scihub + doi)
    if shpage.status_code != 200:
        raise SciHubError('Sci-Hub website unaccessible.')
    if shpage.text.find('article not found') > 0:
        raise ArticleNotFound()
    else:
        shparse = bs4.BeautifulSoup(shpage.text)
        return 'http://' + shparse.body.find(
            'iframe',
            attrs={'id': 'pdf'}
        ).get('src').lstrip('/')

# ------------------

def textbox(master, label, width, textwrap, row, column):
    # make the text box with x and y scrollbar
    frame = tk.Frame(master)
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

# ------------------

def massagebox(master, text):
    mb = tk.Toplevel(master)
    text = tk.Label(mb, text=text)
    text.pack()

# ------------------

def logbox(master, label, text):
    lb = tk.Toplevel(master)
    log = textbox(lb, label, 80, tk.NONE, 0, 0)
    log.insert(tk.END, text)


# ------------------

class Woodpecker(tk.Frame):
    def __init__(self, master):
        self.data = dict()
        super().__init__(master)
        self.pack()
        self.doi_textbox = textbox(
            self,
            'DOI 输入框',
            width=40, textwrap=tk.NONE,
            row=0, column=0
        )
        path_entry_frame = tk.Frame(self)
        path_entry_frame.grid(row=1, column=0)
        self.path_label = tk.Label(
            path_entry_frame, text='保存路径'
        )
        self.path_label.pack(anchor=tk.W)
        self.path_entry = tk.Entry(
            path_entry_frame, width=30
        )
        self.path_entry.pack(anchor=tk.E)
        self.execute_button = tk.Button(
            self,
            text='Woodpecker!',
            command=self.button_execute
        )
        self.execute_button.grid(row=2, column=0)

    def button_execute(self):
        self.data['dois'] = self.doi_textbox.get(
            1.0,
            tk.END
        ).rstrip().lower().split('\n')
        self.data['dir'] = self.path_entry.get()
        self.data['log'] = list()
        for doi in self.data['dois']:
            try:
                try:
                    pdfurl = scihub_pdfurl(doi, goodscihub)
                    pdf = requests.get(pdfurl)
                    pdfname = os.path.join(
                        self.data['dir'],
                        doi.replace('/', '') + '.pdf'
                    )
                    with open(pdfname, 'wb') as f:
                        f.write(pdf.content)
                    self.data['log'].append(
                        ': '.join(
                            [
                                doi,
                                time.strftime('%Y/%m/%d %H:%M'),
                                'OK!'
                            ]
                        )
                    )
                except (ConnectionError, SciHubError):
                    massagebox(
                        self,
                        'Sci-Hub 无法链接，如果确定网络没有问题，那说明 Sci-Hub 暂时被封了，请下载新版本或者联系本软件作者 wolfsonliu@gmail.com。'
                    )
                    break
            except ArticleNotFound:
                self.data['log'].append(
                    ': '.join(
                        [
                            doi,
                            time.strftime('%Y/%m/%d %H:%M'),
                            '[ failed ]'
                        ]
                    )
                )
            finally:
                pass
        with open(
                os.path.join(self.data['dir'], 'woodpecker.log'), 'a'
        ) as f:
            f.write('\n'.join(self.data['log']))
        self.logbox = logbox(
            self, 'Woodpecker 运行状况', '\n'.join(self.data['log'])
        )

# ------------------

if __name__ == '__main__':
    root = tk.Tk()
    root.title('Woodpecker')
    app = Woodpecker(root)
    root.mainloop()

# ------------------
# EOF
# ------------------
