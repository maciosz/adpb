import sys
from Tkinter import *

def format_error():
    toplevel = Toplevel()
    label1 = Label(toplevel, height=8, width=45, text=FERROR_TEXT)
    label1.pack()
    return toplevel
FERROR_TEXT = "Invalid conversion. Change type of input or output file."
