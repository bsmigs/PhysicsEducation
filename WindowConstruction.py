import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import Tkinter as tk


def loadFields(root, fields, unitFields, currState='normal'):
    entries = []
    for key in fields:
        row = tk.Frame(root)
        label = tk.Label(row, width=15, text=key, anchor='w')
        entry = tk.Entry(row)
        entry.insert(tk.END, fields[key])

        row.bind()
        row.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)
        label.pack(side=tk.LEFT)
        entry.pack(side=tk.LEFT, expand=tk.YES, fill=tk.X)
        entry.config(state=currState)

        # set units
        unitsLabel = tk.Label(row, width=8, text=unitFields[key], anchor='w')
        unitsLabel.pack(side=tk.RIGHT)

        # append values to entries list
        entries.append( (key, entry) )

    return entries
