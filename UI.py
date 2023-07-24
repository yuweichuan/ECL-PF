from tkinter import *
from tkinter import filedialog
from tkinter import ttk
from tkinter import messagebox
import threading
import time
from PIL import ImageTk, Image
import subprocess
import webbrowser
import sqlite3
import json

'''############################################## Main page functions ###############################################'''


def on_enter_noncleavable(e):
    """Main page: click noncleavable button"""
    e.widget['background'] = 'white'
    status_label.config(text='The non-cleavable searching module supports the data derived from non-cleavable cross-linkers.')


def on_enter_cleavable(e):
    """Main page: click cleavable button"""
    e.widget['background'] = 'white'
    status_label.config(text='The cleavable searching module supports the data derived from cleavable cross-linkers.')


def on_leave(e):
    """Main page: leave button"""
    e.widget['background'] = 'SystemButtonFace'
    status_label.config(text='Welcome to ECL 3.0!')


def callurl(url):
    """Main page: open url"""
    webbrowser.open_new_tab(url)

'''################################################ Common functions ################################################'''


def browse():
    """Browse data set"""
    global infile
    data_en.configure(state=NORMAL)
    data_en.delete('1.0', END)
    infile = filedialog.askopenfilenames(title = "Select a file",filetypes = (("MZXML File","*.mzXML"),))
    for ele_file in infile:
        data_en.insert(END, ele_file + '\n')
    data_en.configure(state=DISABLED)


def browsefasta():
    """Browse database file (FASTA)"""
    global fasta
    db_en.configure(state=NORMAL)
    db_en.delete(0, END)
    fasta = filedialog.askopenfilename(title = "Select a file",filetypes = (("FASTA","*.fasta"),))
    db_en.insert(0,fasta)
    db_en.configure(state=DISABLED)


def back2normal1():
    """Return to main page"""
    noncleavable['state'] = NORMAL
    cleavable['state'] = NORMAL
    non_win.destroy()
    root.deiconify()


def handle_focus_in(key):
    """Click and show"""
    key.delete(0, END)
    key.config(fg='black')


def addfixedmod():
    """add fix modification"""
    if mod_en.curselection():
        item = mod_en.get(ANCHOR)
        fixed_en.insert(END, item)


def addvarmod():
    """add variable modification"""
    if mod_en.curselection():
        item = mod_en.get(ANCHOR)
        var_en.insert(END, item)


def remfixedmod():
    """remove fix modification"""
    if fixed_en.curselection():
        index = fixed_en.curselection()[0]
        fixed_en.delete(index)


def remvarmod():
    """remove variable modification"""
    if var_en.curselection():
        index = var_en.curselection()[0]
        var_en.delete(index)


def stop():
    """Stop running"""
    global flag1
    response = messagebox.askyesno('Stop task', 'Are you sure to stop the running task?')
    if response == 1:
        print_en.insert(END, '\n{}Task terminated!{}\n'.format('*'*47, '*'*47))
        flag1 = True
        p1.kill()
        # print_en.configure(state=DISABLED)
        stop_btn['state']=DISABLED
        messagebox.showinfo('Task', 'Task terminated!')


def on_enter_data_btn(e):
    """Enter dataset browse button"""
    e.widget['background'] = 'white'
    non_status.config(text='Browse data file(s) in .mzXML format (centroid MS1+MS2 spectra)')


def on_enter_db_btn(e):
    """Enter database browse button"""
    e.widget['background'] = 'white'
    non_status.config(text='Browse protein database in .FASTA format (decoy sequence excluded)')


def on_enter_edit_btn(e):
    """Enter edit linker button"""
    e.widget['background'] = 'white'
    non_status.config(text='Edit the current selected cross linker in the ECL database')


def on_enter_add_btn(e):
    """Enter add linker button"""
    e.widget['background'] = 'white'
    non_status.config(text='Add new cross linker in the ECL database')


def on_enter_del_btn(e):
    """Enter delete linker button"""
    e.widget['background'] = 'white'
    non_status.config(text='Delete the current selected cross linker in the ECL database')


def on_enter_ref_btn(e):
    """Enter linker refresh button"""
    e.widget['background'] = 'white'
    non_status.config(text='Update the ECL cross linker database each time after the add, edit or delete operation')


def on_enter_mod_edit_btn(e):
    """Enter edit modification button"""
    e.widget['background'] = 'white'
    non_status.config(text='Select a modification and edit it in the ECL database')


def on_enter_mod_add_btn(e):
    """Enter edit modification button"""
    e.widget['background'] = 'white'
    non_status.config(text='Add a new modification in the ECL database')


def on_enter_mod_del_btn(e):
    """Enter delete modification button"""
    e.widget['background'] = 'white'
    non_status.config(text='Select a modification and delete it in the ECL database')


def on_enter_mod_ref_btn(e):
    """Enter refresh modification button"""
    e.widget['background'] = 'white'
    non_status.config(text='Update the ECL modification database each time after the add, edit or delete operation')


def on_enter_add_fixed_btn(e):
    """Enter fix modification button"""
    e.widget['background'] = 'white'
    non_status.config(text='Select a mod and add it as a fixed modification in the task setting')


def on_enter_add_var_btn(e):
    """Enter variable modification button"""
    e.widget['background'] = 'white'
    non_status.config(text='Select a mod and add it as a variable modification in the task setting')


def on_enter_rem_fixed_btn(e):
    """Enter remove fix modification button"""
    e.widget['background'] = 'white'
    non_status.config(text='Select a mod and remove it from the fixed modification')


def on_enter_rem_var_btn(e):
    """Enter remove variable modification button"""
    e.widget['background'] = 'white'
    non_status.config(text='Select a mod and remove it from the variable modification')


def on_enter_start_btn(e):
    """Enter start button"""
    e.widget['background'] = 'white'
    non_status.config(text='Start the task running')


def on_enter_stop_btn(e):
    """Enter stop button"""
    e.widget['background'] = 'white'
    non_status.config(text='Stop the task running')


def on_enter_back_btn(e):
    """Enter back button"""
    e.widget['background'] = 'white'
    non_status.config(text='Back to the main page')


def on_enter_linker_name_en(e):
    """Enter linker name button"""
    add_status.config(text='Input cross linker name')


def on_enter_linker_site_en(e):
    """Enter linker site button"""
    add_status.config(text='Input one amino acid name in capital letter')


def on_enter_linker_mass_en(e):
    """Enter mono mass button"""
    add_status.config(text='Input cross linker mono (residual) mass')


def on_enter_mod_name_en(e):
    """Enter modification name"""
    add_status.config(text='Input mods. name with only lowercase letters allowed')


def on_enter_mod_site_en(e):
    """Enter modified amino acid"""
    add_status.config(text='Input one amino acid name in capital letter')


def on_enter_mod_mass_en(e):
    """Enter modification mass"""
    add_status.config(text='Input modification mono mass')
'''################################################ Cleavable frame #################################################'''


def stop_cleave():
    """Cleavable module: stop function"""
    global flag2
    response = messagebox.askyesno('Stop task', 'Are you sure to stop the running task?')
    if response == 1:

        flag2 = True
        if p1.poll() is None:
            p1.kill()
        elif p2.poll() is None:
            p2.kill()
        elif align_var.get() == 1 and p3.poll() is None:
            p3.kill()
        elif precursor_var.get() == 1 and p4.poll() is None:
            p4.kill()
        elif p5.poll() is None:
            p5.kill()
        elif p6.poll() is None:
            p6.kill()

        # print_en.configure(state=DISABLED)
        print_en.insert(END, '\n{}Task terminated!{}\n'.format('*' * 47, '*' * 47))
        stop_btn['state']=DISABLED
        messagebox.showinfo('Task', 'Task terminated!')


def cleave_params_print():
    """Cleavable module: print parameters and generate ECLPF_conf"""
    conf = dict()
    print_en.configure(state=NORMAL)
    print_en.delete(1.0, END)
    print_en.insert(END, "{}Parameters{}\n\n".format('*'*50, '*'*50))
    print_en.insert(END, "Data :              ")
    print_en.insert(END, infile[0] + '\n')
    for ele in infile[1:]:
        print_en.insert(END,'                    ' + ele + '\n')
    conf['data_path'] = infile
    print_en.insert(END, 'FASTA file :        ' + db_en.get() + '\n')
    conf['fasta_path'] = db_en.get()
    print_en.insert(END, "Output:             ")
    print_en.insert(END, infile[0].replace('.mzXML', '_final.csv') + '\n')
    for ele in infile[1:]:
        print_en.insert(END, '                    ' + ele.replace('.mzXML', '_final.csv') + '\n')
    print_en.insert(END, 'Threads :           ' + core_en.get() + '\n')
    conf['thread'] = int(core_en.get())
    print_en.insert(END, 'Peptide length :    ' + '[' + lower_en.get() + ', ' + higher_en.get() + ']' + '\n')
    conf['max_length'] = int(higher_en.get())
    conf['min_length'] = int(lower_en.get())
    print_en.insert(END, 'Digestion :         ' + parse_en.get() + '\n')
    conf['parse_rule'] = parse_en.get()
    print_en.insert(END, 'Miss cleavage :     ' + miss_en.get() + '\n')
    conf['miss_cleavage'] = int(miss_en.get())
    print_en.insert(END, 'Ms1 tolerance :     ' + ms1tol_en.get() + ' ppm\n')
    conf['ms1_tol'] = round(int(ms1tol_en.get()) * 0.000001, 8)
    print_en.insert(END, 'Ms2 tolerance :     ' + ms2tol_en.get() + ' ppm\n')
    conf['ms2_tol'] = round(int(ms2tol_en.get()) * 0.000001, 8)
    print_en.insert(END, 'Cross linker :      ' + cleavelinker_en.get() + '\n')
    print_en.insert(END, 'Link site :         ' + linkerDict[cleavelinker_en.get()][1] + '\n')
    conf['link_site'] = [xls for xls in linkerDict[cleavelinker_en.get()][1]]
    print_en.insert(END, 'Linker mass :       ' + str(linkerDict[cleavelinker_en.get()][2]) + '\n')
    conf['xl_mass'] = float(linkerDict[cleavelinker_en.get()][2])
    print_en.insert(END, 'Short res. mass :   ' + str(linkerDict[cleavelinker_en.get()][3]) + '\n')
    conf['m_short'] = float(linkerDict[cleavelinker_en.get()][3])
    print_en.insert(END, 'Long res. mass :    ' + str(linkerDict[cleavelinker_en.get()][4]) + '\n')
    conf['m_long'] = float(linkerDict[cleavelinker_en.get()][4])

    conf['fix_mod'] = dict()
    if fixed_en.get(0, END):
        print_en.insert(END, 'Fixed mods. :       ' + fixed_en.get(0, END)[0] + ' ' + str(modDict[fixed_en.get(0, END)[0]][2]) + '\n')
        conf['fix_mod'][modDict[fixed_en.get(0, END)[0]][0]] = [modDict[fixed_en.get(0, END)[0]][2],
                                                                modDict[fixed_en.get(0, END)[0]][1].split('&')]
        for ele in fixed_en.get(0, END)[1:]:
            print_en.insert(END, '                    ' + ele + ' ' + str(modDict[ele][2]) + '\n')
            conf['fix_mod'][modDict[ele][0]] = [modDict[ele][2], modDict[ele][1].split('&')]

    conf['var_mod'] = dict()
    if var_en.get(0, END):
        print_en.insert(END, 'Var mods. :         ' + var_en.get(0, END)[0] + ' ' + str(modDict[var_en.get(0, END)[0]][2]) + '\n')
        conf['var_mod'][modDict[var_en.get(0, END)[0]][0]] = [modDict[var_en.get(0, END)[0]][2],
                                                              modDict[var_en.get(0, END)[0]][1].split('&')]
        for ele in var_en.get(0, END)[1:]:
            print_en.insert(END, '                    ' + ele + ' ' + str(modDict[ele][2]) + '\n')
            conf['var_mod'][modDict[ele][0]] = [modDict[ele][2], modDict[ele][1].split('&')]
    print_en.insert(END, 'Max mods./peptide : ' + max_mod_en.get() + '\n')
    conf['num_max_mod'] = int(max_mod_en.get())
    print_en.insert(END, 'MS2 act. type :     ' + activation_en.get() + '\n')
    if activation_en.get() == '[HCD]' or activation_en.get() == '[HCD, ETD]':
        conf['activation_type'] = ['HCD', 'ETD']
    else:
        conf['activation_type'] = ['CID', 'ETD']
    print_en.insert(END, 'Prec. mass refine : ' + ('True' if precursor_var.get() else 'False') + '\n')
    print_en.insert(END, 'Local alignment :   ' + ('True' if align_var.get() else 'False') + '\n')
    print_en.insert(END, 'FDR threshold :     ' + fdr_en.get() + ' %\n\n')
    print_en.insert(END, "{}Program running{}\n\n".format('*' *47, '*' * 48))

    with open('ECLPF_conf', 'w') as f:
        json.dump(conf, f, indent=4, separators=(',', ':'))


def sub_cleave():
    """Cleavable module: start function"""
    global precursor_var, align_var, fdr_en, flag2, p1, p2, p3, p4, p5, p6
    flag2 = False
    print_en.insert(END, '>> Begin at {}\n'.format(time.ctime()))
    print_en.insert(END, '>> Generating protein sequence database...\n')
    p1 = subprocess.Popen(['python', 'database.py'], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          stdin=subprocess.PIPE)

    p1.wait()
    if flag2:
        return True
    time.sleep(1)

    print_en.insert(END, '>> Preparing data...\n')
    p2 = subprocess.Popen(['python', 'spectra_separation.py'], shell=True, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, stdin=subprocess.PIPE)

    p2.wait()
    if flag2:
        return True
    time.sleep(1)

    if align_var.get() == 1:
        print_en.insert(END, '>> Local alignment processing...\n')
        p3 = subprocess.Popen(['python', 'local_alignment.py'], shell=True, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, stdin=subprocess.PIPE)
        p3.wait()
        if flag2:
            return True
        time.sleep(1)

    if precursor_var.get() == 1:
        print_en.insert(END, '>> Precursor mass refining...\n')
        p4 = subprocess.Popen(['python', 'precursor_refinement.py'], shell=True, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, stdin=subprocess.PIPE)
        while p4.poll() is None:
            for line4 in iter(p4.stdout.readline, b""):
                print_en.insert(END, line4.decode("utf-8"))
                print_en.see(END)
        p4.wait()
        if flag2:
            return True
        time.sleep(1)

    print_en.insert(END, '>> Searching ...\n')
    p5 = subprocess.Popen(['python', 'ECL_PF.py'], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          stdin=subprocess.PIPE)

    p5.wait()
    if flag2:
        return True
    time.sleep(1)

    print_en.insert(END, '>> Filtering ...\n')
    with open('ECLPF_conf', 'r') as conf_file:
        paths = json.load(conf_file)['data_path']

    for ele_path in paths:
        inpath = ele_path.split('/')[-1].replace('.mzXML', '.csv')
        outpath = ele_path.replace('.mzXML', '_final.csv')
        command_f = 'python cleavable_splitctrl.py {} {} {}'.format(inpath, outpath, round(float(fdr_en.get())*0.01, 2))
        p6 = subprocess.Popen(command_f, shell=False, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, stdin=subprocess.PIPE)
        p6.wait()
        if flag2:
            return True

    print_en.insert(END, '>> End at {}\n\n'.format(time.ctime()))
    print_en.insert(END, "{}Program finished{}\n".format('*' * 47, '*' * 47))

    print_en.configure(state=DISABLED)
    stop_btn['state'] = DISABLED
    messagebox.showinfo('Task', 'Task finished!')


def run_cleave():
    """Cleavable module: check input and run"""
    global infile
    """check input params first and then call funtions to run program"""
    if not infile:
        messagebox.showwarning('Invalid parameter settings', 'No input data file(s)')
        return True
    for ele in infile:
        if ' ' in ele:
            messagebox.showwarning('Invalid parameter settings', 'Space is contained in the path')
            return True
    infile = list(infile)

    if not db_en.get():
        messagebox.showwarning('Invalid parameter settings', 'No input FASTA file')
        return True

    try:
        if float(lower_en.get()) > float(higher_en.get()):
            messagebox.showwarning('Invalid parameter settings', 'Invalid peptide length restriction')
            return True
    except ValueError:
        messagebox.showwarning('Invalid parameter settings', 'Invalid peptide length')
        return True

    try:
        if float(ms1tol_en.get()) <=0 or float(ms2tol_en.get()) <= 0:
            messagebox.showwarning('Invalid parameter settings', 'Mass tolerance is smaller than 0')
            return True
    except ValueError:
        messagebox.showwarning('Invalid parameter settings', 'Invalid mass tolerance')
        return True

    try:
        if float(fdr_en.get()) <= 0:
            messagebox.showwarning('Invalid parameter settings', 'FDR is less than 0')
            return True
    except ValueError:
        messagebox.showwarning('Invalid parameter settings', 'Invalid FDR setting')
        return True

    stop_btn['state'] = NORMAL
    print_en.configure(state=NORMAL)
    cleave_params_print()

    threading.Thread(target=sub_cleave).start()


def cleav():
    """Cleavable module: frame"""
    global non_win
    global img_back, img_start, img_stop
    global data_en, db_en, core_en, infile
    global lower_en, higher_en, ms1tol_en, ms2tol_en, parse_en, miss_en, fdr_en, precursor_var, align_var, activation_en
    global cleavelinker_en, linkerDict
    global fixed_en, var_en, mod_en, modDict, max_mod_en
    global print_en, stop_btn
    global non_status
    noncleavable['state'] = DISABLED
    cleavable['state'] = DISABLED
    root.withdraw()

    non_win = Toplevel()
    non_win.title('Cleavable searching module')
    non_win.iconbitmap('ECL.ico')
    non_win_w = 900
    non_win_h = 700
    non_win_width = non_win.winfo_screenwidth()
    non_win_height = non_win.winfo_screenheight()
    non_win_x = int(non_win_width / 2 - non_win_w / 2)
    non_win_y = int(non_win_height / 2 - non_win_h / 2)
    non_win.geometry('%dx%d+%d+%d' % (non_win_w, non_win_h, non_win_x, non_win_y))
    non_win.protocol("WM_DELETE_WINDOW", back2normal1)
    '''Input frame'''
    input_frame = LabelFrame(non_win, text='Input')
    input_frame.grid(row=0, column=0, rowspan=3, columnspan=2, padx=5, pady=(10, 5), sticky='nsew')
    input_frame.grid_propagate(False)

    data_label = Label(input_frame, text="Data file(s)")
    data_label.grid(row=0, column=0, rowspan=4, sticky="nw")

    data_en_scroll = Scrollbar(input_frame, orient='horizontal')
    data_en_scroll.grid(row=3, column=1, sticky="nsew")
    data_en = Text(input_frame, state=DISABLED, width=10, height=5, wrap=NONE, font=('Arial', 8),
                   xscrollcommand=data_en_scroll.set)
    data_en.grid(row=0, column=1, rowspan=3, sticky="nsew")
    data_en_scroll.config(command=data_en.xview)

    infile = tuple()
    data_btn = Button(input_frame, text="Browse", command=browse)
    data_btn.grid(row=0, column=2, rowspan=4, padx=5, sticky="new")
    data_btn.bind("<Enter>", on_enter_data_btn)
    data_btn.bind("<Leave>", on_cleave_leave)

    db_label = Label(input_frame, text="Database")
    db_label.grid(row=4, column=0, sticky="nsew")
    db_en = Entry(input_frame, font=('Arial', 8), state=DISABLED)
    db_en.grid(row=4, column=1, sticky="nsew")
    db_btn = Button(input_frame, text="Browse", command=browsefasta)
    db_btn.grid(row=4, column=2, padx=5, sticky="nsew")
    db_btn.bind("<Enter>", on_enter_db_btn)
    db_btn.bind("<Leave>", on_cleave_leave)

    core_label = Label(input_frame, text="# threads")
    core_label.grid(row=5, column=0, sticky="nsew")
    core_en = ttk.Combobox(input_frame, width=4, state="readonly")
    core_en['values'] = (1, 2, 4, 8, 16)
    core_en.current(2)
    core_en.grid(row=5, column=1, pady=5, sticky="nsw")

    input_frame.grid_columnconfigure(1, weight=6)
    input_frame.grid_columnconfigure(2, weight=1)
    input_frame.grid_rowconfigure(0, weight=1)
    input_frame.grid_rowconfigure(1, weight=1)
    input_frame.grid_rowconfigure(2, weight=1)

    '''searching frame'''
    search_frame = LabelFrame(non_win, text='Parameters')
    search_frame.grid(row=3, column=0, rowspan=3, columnspan=2, padx=5, pady=(5, 5), sticky='nsew')
    search_frame.grid_propagate(False)

    subframe1 = Frame(search_frame)
    subframe1.grid(row=0, column=0, pady=(10, 5), sticky='nsew')
    subframe1.grid_propagate(False)
    mass_label = Label(subframe1, text="<= peptide length <=")
    mass_label.grid(row=0, column=1, padx=10, sticky='new')
    lower_en = Entry(subframe1, width=2)
    lower_en.insert(0, 6)
    lower_en.bind("<FocusIn>", lambda event: handle_focus_in(key=lower_en))
    lower_en.grid(row=0, column=0, padx=(10, 0), sticky='new')

    higher_en = Entry(subframe1, width=3)
    higher_en.insert(0, 30)
    higher_en.bind("<FocusIn>", lambda event: handle_focus_in(key=higher_en))
    higher_en.grid(row=0, column=2, sticky='new')

    parse_label = Label(subframe1, text="Digestion rule :")
    parse_en = ttk.Combobox(subframe1, width=12, state="readonly")
    parse_en['values'] = ('trypsin', 'arg-c', 'asp-n', 'pepsin ph1.3', 'pepsin ph2.0')
    parse_en.current(0)
    parse_label.grid(row=0, column=3, padx=(40, 5), sticky='new')
    parse_en.grid(row=0, column=4, padx=(0, 10), sticky='new')
    subframe1.grid_columnconfigure(0, weight=1)
    subframe1.grid_columnconfigure(1, weight=1)
    subframe1.grid_columnconfigure(2, weight=1)
    subframe1.grid_columnconfigure(3, weight=1)
    subframe1.grid_columnconfigure(4, weight=1)

    subframe2 = Frame(search_frame)
    subframe2.grid(row=1, column=0, pady=5, sticky='nsew')
    subframe2.grid_propagate(False)
    miss_label = Label(subframe2, text="Miss cleavage :")
    miss_en = ttk.Combobox(subframe2, width=2, state="readonly")
    miss_en['values'] = (0, 1, 2, 3)
    miss_en.current(2)
    miss_label.grid(row=0, column=0, padx=(5, 5), sticky='new')
    miss_en.grid(row=0, column=1, sticky='new')

    ms1tol_label = Label(subframe2, text="MS1 tolerance :")
    ms1tol_en = Entry(subframe2, width=4)
    ms1tol_en.insert(0, 10)
    ms1unit_label = Label(subframe2, text="ppm")
    ms1tol_label.grid(row=0, column=2, padx=(10, 5), sticky='new')
    ms1tol_en.grid(row=0, column=3, sticky='new')
    ms1tol_en.bind("<FocusIn>", lambda event: handle_focus_in(key=ms1tol_en))
    ms1unit_label.grid(row=0, column=4, sticky='new')

    ms2tol_label = Label(subframe2, text="MS2 tolerance :")
    ms2tol_en = Entry(subframe2, width=4)
    ms2tol_en.insert(0, 20)
    ms2unit_label = Label(subframe2, text="ppm")
    ms2tol_label.grid(row=0, column=5, padx=(10, 5), sticky='new')
    ms2tol_en.grid(row=0, column=6, sticky='new')
    ms2tol_en.bind("<FocusIn>", lambda event: handle_focus_in(key=ms2tol_en))
    ms2unit_label.grid(row=0, column=7, padx=(0, 10), sticky='new')

    subframe2.grid_columnconfigure(0, weight=1)
    subframe2.grid_columnconfigure(1, weight=1)
    subframe2.grid_columnconfigure(2, weight=1)
    subframe2.grid_columnconfigure(3, weight=1)
    subframe2.grid_columnconfigure(4, weight=1)
    subframe2.grid_columnconfigure(5, weight=1)
    subframe2.grid_columnconfigure(6, weight=1)
    subframe2.grid_columnconfigure(7, weight=1)

    subframe3 = Frame(search_frame)
    subframe3.grid(row=2, column=0, pady=5, sticky='nsew')
    subframe3.grid_propagate(False)

    activation_label = Label(subframe3, text="MS2 Activation :")
    activation_en = ttk.Combobox(subframe3, width=10, state="readonly")
    activation_en['values'] = ('[HCD, ETD]', '[CID, ETD]', '[HCD]', '[CID]')
    activation_en.current(0)
    activation_label.grid(row=0, column=2, padx=(5, 5), sticky='new')
    activation_en.grid(row=0, column=3, padx=(0, 5), sticky='new')

    precursor_var = IntVar()
    precursor_en = Checkbutton(subframe3, text='Precursor refinement', variable=precursor_var)
    precursor_en.select()
    precursor_en.grid(row=0, column=0, sticky='new')

    align_var = IntVar()
    align_en = Checkbutton(subframe3, text='Local alignment', variable=align_var)
    align_en.grid(row=0, column=1, sticky='new')

    subframe3.grid_columnconfigure(0, weight=1)
    subframe3.grid_columnconfigure(1, weight=1)
    subframe3.grid_columnconfigure(2, weight=1)
    subframe3.grid_columnconfigure(3, weight=1)

    subframe4 = Frame(search_frame)
    subframe4.grid(row=3, column=0, pady=5, sticky='nsew')
    fdr_label = Label(subframe4, text="FDR setting :")
    fdr_en = Entry(subframe4, width=2)
    fdr_en.insert(0, 1)
    fdr_en.bind("<FocusIn>", lambda event: handle_focus_in(key=fdr_en))

    fdrunit_label = Label(subframe4, text="%")
    fdr_label.grid(row=0, column=0, padx=5)
    fdr_en.grid(row=0, column=1)
    fdrunit_label.grid(row=0, column=2)

    search_frame.grid_columnconfigure(0, weight=1)
    search_frame.grid_rowconfigure(0, weight=1)
    search_frame.grid_rowconfigure(1, weight=1)
    search_frame.grid_rowconfigure(2, weight=1)
    search_frame.grid_rowconfigure(3, weight=1)

    '''linker frame'''
    linkerDict = {}
    linker_frame = LabelFrame(non_win, text='Cross-linker')
    linker_frame.grid(row=0, column=2, rowspan=2, columnspan=2, padx=5, pady=(10, 5), sticky='nsew')
    linker_frame.grid_propagate(False)
    cleavelinker_label = Label(linker_frame, text='Select cleavable cross linker :')
    cleavelinker_en = ttk.Combobox(linker_frame, width=4, state="readonly")
    refresh_cleavelinker()
    cleavelinker_label.grid(row=0, column=0, columnspan=3, padx=5, pady=5, sticky='nsew')
    cleavelinker_en.grid(row=0, column=3, padx=20, pady=10, sticky='nsew')
    db_frame = LabelFrame(linker_frame, text='Cleavable cross linker database')
    db_frame.grid(row=1, column=0, columnspan=3, padx=5, pady=(10, 5), sticky='nsew')
    db_frame.grid_propagate(False)
    add_btn = Button(db_frame, text='Add', width=8, command=add_cleavelinker)
    edit_btn = Button(db_frame, text='Edit', width=8, command=edit_cleavelinker)
    del_btn = Button(db_frame, text='Delete', width=8, command=delete_cleavelinker)
    ref_btn = Button(linker_frame, text='Update database', width=10, command=refresh_cleavelinker)
    add_btn.grid(row=0, column=0, padx=10, pady=5, sticky='nsew')
    add_btn.bind("<Enter>", on_enter_add_btn)
    add_btn.bind("<Leave>", on_cleave_leave)
    edit_btn.grid(row=0, column=1, padx=10, pady=5, sticky='nsew')
    edit_btn.bind("<Enter>", on_enter_edit_btn)
    edit_btn.bind("<Leave>", on_cleave_leave)
    del_btn.grid(row=0, column=2, padx=10, pady=5, sticky='nsew')
    del_btn.bind("<Enter>", on_enter_del_btn)
    del_btn.bind("<Leave>", on_cleave_leave)
    ref_btn.grid(row=1, column=3, padx=20, pady=(20, 10), sticky='nsew')
    ref_btn.bind("<Enter>", on_enter_ref_btn)
    ref_btn.bind("<Leave>", on_cleave_leave)

    db_frame.grid_columnconfigure(0, weight=1)
    db_frame.grid_columnconfigure(1, weight=1)
    db_frame.grid_columnconfigure(2, weight=1)
    db_frame.grid_rowconfigure(0, weight=1)

    linker_frame.grid_rowconfigure(0, weight=1)
    linker_frame.grid_rowconfigure(1, weight=1)
    linker_frame.grid_columnconfigure(0, weight=1)
    linker_frame.grid_columnconfigure(1, weight=1)
    linker_frame.grid_columnconfigure(2, weight=1)
    linker_frame.grid_columnconfigure(3, weight=1)

    '''mod frame'''
    mod_frame = LabelFrame(non_win, text='Modifications')
    mod_frame.grid(row=2, column=2, rowspan=4, columnspan=2, padx=5, pady=(5, 5), sticky='nsew')
    mod_frame.grid_propagate(False)
    mod_db_frame = LabelFrame(mod_frame, text='Mods. database')
    mod_db_frame.grid(row=0, column=3, rowspan=3, padx=5, pady=(10, 5), sticky='nsew')
    mod_db_frame.grid_propagate(False)

    fixed_label = Label(mod_frame, text="fixed mods.")
    fixed_en = Listbox(mod_frame, height=7, width=15)
    var_label = Label(mod_frame, text="variable mods.")
    var_en = Listbox(mod_frame, height=7, width=15)
    fixed_label.grid(row=0, column=0, pady=5, sticky='nsew')
    fixed_en.grid(row=1, column=0, padx=5, pady=5, sticky='nsew')
    var_label.grid(row=2, column=0, pady=5, sticky='nsew')
    var_en.grid(row=3, column=0, padx=5, pady=5, sticky='nsew')

    add_fixed_btn = Button(mod_frame, text="<", command=addfixedmod, anchor='n')
    rem_fixed_btn = Button(mod_frame, text=">", command=remfixedmod, anchor='s')
    add_var_btn = Button(mod_frame, text="<", command=addvarmod, anchor='n')
    rem_var_btn = Button(mod_frame, text=">", command=remvarmod, anchor='s')
    add_fixed_btn.grid(row=1, column=1, pady=15, sticky='n')
    add_fixed_btn.bind("<Enter>", on_enter_add_fixed_btn)
    add_fixed_btn.bind("<Leave>", on_cleave_leave)
    rem_fixed_btn.grid(row=1, column=1, pady=15, sticky='s')
    rem_fixed_btn.bind("<Enter>", on_enter_rem_fixed_btn)
    rem_fixed_btn.bind("<Leave>", on_cleave_leave)
    add_var_btn.grid(row=3, column=1, pady=15, sticky='n')
    add_var_btn.bind("<Enter>", on_enter_add_var_btn)
    add_var_btn.bind("<Leave>", on_cleave_leave)
    rem_var_btn.grid(row=3, column=1, pady=15, sticky='s')
    rem_var_btn.bind("<Enter>", on_enter_rem_var_btn)
    rem_var_btn.bind("<Leave>", on_cleave_leave)

    max_mod_label = Label(mod_frame, text="Max mods./pep")
    max_mod_en = ttk.Combobox(mod_frame, width=2, state="readonly")
    max_mod_en['values'] = (0, 1, 2, 3, 4)
    max_mod_en.current(2)
    max_mod_label.grid(row=0, column=2, padx=5, sticky='w')
    max_mod_en.grid(row=0, column=2, padx=5, sticky='e')
    modDict = {}
    mod_en = Listbox(mod_frame, height=18, width=15)
    refresh_mod_cleavelinker()
    mod_en.grid(row=1, column=2, rowspan=3, padx=5, pady=5, sticky='nsew')

    mod_add_btn = Button(mod_db_frame, text='Add', width=8, command=add_mod_cleavelinker)
    mod_edit_btn = Button(mod_db_frame, text='Edit', width=8, command=edit_mod_cleavelinker)
    mod_del_btn = Button(mod_db_frame, text='Delete', width=8, command=delete_mod_cleavelinker)
    mod_ref_btn = Button(mod_frame, text='Update mods.', width=10, command=refresh_mod_cleavelinker)
    mod_add_btn.grid(row=0, column=0, padx=5, pady=5, sticky='nsew')
    mod_add_btn.bind("<Enter>", on_enter_mod_add_btn)
    mod_add_btn.bind("<Leave>", on_cleave_leave)
    mod_edit_btn.grid(row=1, column=0, padx=5, pady=5, sticky='nsew')
    mod_edit_btn.bind("<Enter>", on_enter_mod_edit_btn)
    mod_edit_btn.bind("<Leave>", on_cleave_leave)
    mod_del_btn.grid(row=2, column=0, padx=5, pady=5, sticky='nsew')
    mod_del_btn.bind("<Enter>", on_enter_mod_del_btn)
    mod_del_btn.bind("<Leave>", on_cleave_leave)
    mod_ref_btn.grid(row=3, column=3, padx=5, pady=(25, 25), sticky='nsew')
    mod_ref_btn.bind("<Enter>", on_enter_mod_ref_btn)
    mod_ref_btn.bind("<Leave>", on_cleave_leave)

    mod_db_frame.grid_columnconfigure(0, weight=1)
    mod_db_frame.grid_rowconfigure(0, weight=1)
    mod_db_frame.grid_rowconfigure(1, weight=1)
    mod_db_frame.grid_rowconfigure(2, weight=1)

    mod_frame.grid_rowconfigure(0, weight=1)
    mod_frame.grid_rowconfigure(1, weight=3)
    mod_frame.grid_rowconfigure(2, weight=1)
    mod_frame.grid_rowconfigure(3, weight=3)
    mod_frame.grid_columnconfigure(0, weight=2)
    mod_frame.grid_columnconfigure(1, weight=1)
    mod_frame.grid_columnconfigure(2, weight=3)
    mod_frame.grid_columnconfigure(3, weight=2)

    '''back, stop, start'''

    img_back = ImageTk.PhotoImage(Image.open(r"arrow.png").resize((25, 20)))
    back_btn = Button(non_win, text="Back", image=img_back, borderwidth=2, compound=LEFT, command=back2normal1)
    back_btn.grid(row=6, column=0, padx=5, ipadx=10, ipady=5, sticky='w')
    back_btn.bind("<Enter>", on_enter_back_btn)
    back_btn.bind("<Leave>", on_cleave_leave)
    img_stop = ImageTk.PhotoImage(Image.open(r"stop.png").resize((20, 20)))
    stop_btn = Button(non_win, text='Stop', image=img_stop, borderwidth=2, compound=LEFT, state=DISABLED, command=stop_cleave)
    stop_btn.grid(row=6, column=3, padx=5, ipadx=10, ipady=5, sticky='w')
    stop_btn.bind("<Enter>", on_enter_stop_btn)
    stop_btn.bind("<Leave>", on_cleave_leave)
    img_start = ImageTk.PhotoImage(Image.open(r"start.png").resize((25, 20)))
    start_btn = Button(non_win, text='Start', image=img_start, borderwidth=2, compound=LEFT, command=run_cleave)
    start_btn.grid(row=6, column=3, padx=5, ipadx=10, ipady=5, sticky='e')
    start_btn.bind("<Enter>", on_enter_start_btn)
    start_btn.bind("<Leave>", on_cleave_leave)

    '''print frame'''
    print_frame = Frame(non_win)
    print_frame.grid(row=7, column=0, rowspan=3, columnspan=4, padx=5, pady=(5, 5), sticky='nsew')
    print_frame.grid_propagate(False)
    print_en = Text(print_frame, bg="white", relief="sunken", wrap=NONE, state=DISABLED)
    print_en.grid(row=0, column=0, sticky='nsew')

    print_frame.grid_columnconfigure(0, weight=1)
    print_frame.grid_rowconfigure(0, weight=1)
    '''status'''
    non_status = Label(non_win, text='Cleavable cross-linker searching module', bd=1, relief=SUNKEN, anchor=E)
    non_status.grid(row=10, column=0, columnspan=4, ipady=4, sticky='sew')

    non_win.grid_columnconfigure(0, weight=1)
    non_win.grid_columnconfigure(1, weight=1)
    non_win.grid_columnconfigure(2, weight=1)
    non_win.grid_columnconfigure(3, weight=1)

    non_win.grid_rowconfigure(0, weight=1)
    non_win.grid_rowconfigure(1, weight=1)
    non_win.grid_rowconfigure(2, weight=1)
    non_win.grid_rowconfigure(3, weight=1)
    non_win.grid_rowconfigure(4, weight=1)
    non_win.grid_rowconfigure(5, weight=1)
    # non_win.grid_rowconfigure(6, weight=1)
    non_win.grid_rowconfigure(7, weight=1)
    non_win.grid_rowconfigure(8, weight=1)
    non_win.grid_rowconfigure(9, weight=1)
    # non_win.grid_rowconfigure(10, weight=1)


def save_cleave_add():
    """Cleavable module: save added linker"""
    if not linker_en.get():
        messagebox.showwarning('Invalid input', 'No linker name')
        return True
    if len(site_en.get()) != 1 or not site_en.get().isupper():
        messagebox.showwarning('Invalid input', 'Invalid link site')
        return True
    try:
        float(mass_en.get())
    except ValueError:
        messagebox.showwarning('Invalid input', 'Invalid linker mass')
        return True

    try:
        float(short_en.get())
    except ValueError:
        messagebox.showwarning('Invalid input', 'Invalid short residue mass')
        return True

    try:
        float(long_en.get())
    except ValueError:
        messagebox.showwarning('Invalid input', 'Invalid long residue mass')
        return True

    conn = sqlite3.connect('paramsDB.db')
    c = conn.cursor()
    c.execute("INSERT INTO cleavelinker VALUES (:linker, :site, :mass, :m_short, :m_long)",
              {
                  'linker': linker_en.get(),
                  'site': site_en.get(),
                  'mass': mass_en.get(),
                  'm_short': short_en.get(),
                  'm_long': long_en.get()
              })
    conn.commit()
    conn.close()
    add_win.destroy()


def save_cleave_edit():
    """Cleavable module: save edited linker"""
    if not linker_en.get():
        messagebox.showwarning('Invalid input', 'No linker name')
        return True
    if len(site_en.get()) != 1 or not site_en.get().isupper():
        messagebox.showwarning('Invalid input', 'Invalid link site')
        return True
    try:
        float(mass_en.get())
    except ValueError:
        messagebox.showwarning('Invalid input', 'Invalid linker mass')
        return True

    try:
        float(short_en.get())
    except ValueError:
        messagebox.showwarning('Invalid input', 'Invalid short residue mass')
        return True

    try:
        float(long_en.get())
    except ValueError:
        messagebox.showwarning('Invalid input', 'Invalid long residue mass')
        return True

    conn = sqlite3.connect('paramsDB.db')
    c = conn.cursor()
    c.execute("""UPDATE cleavelinker SET
        linker = :linker,
        site = :site,
        mass = :mass,
        m_short = :m_short,
        m_long = :m_long
        WHERE oid = :oid""",
              {'linker': linker_en.get(),
               'site': site_en.get(),
               'mass': mass_en.get(),
               'm_short': short_en.get(),
               'm_long': long_en.get(),
               'oid': linkerDict[cleavelinker_en.get()][5]})
    conn.commit()
    conn.close()
    add_win.destroy()


def add_cleavelinker():
    """Cleavable module: add linker"""
    global linker_en, site_en, mass_en, short_en, long_en, add_win, add_status
    add_win = Toplevel()
    add_win.title('Add cleavable linker')
    add_win.iconbitmap('ECL.ico')
    add_win_w = 300
    add_win_h = 200
    add_win_width = add_win.winfo_screenwidth()
    add_win_height = add_win.winfo_screenheight()
    add_win_x = int(add_win_width / 2 - add_win_w / 2)
    add_win_y = int(add_win_height / 2 - add_win_h / 2)
    add_win.geometry('%dx%d+%d+%d' % (add_win_w, add_win_h, add_win_x, add_win_y))

    linker_en = Entry(add_win, width=20)
    linker_label = Label(add_win, text='Cross linker name :')
    site_en = Entry(add_win, width=20)
    site_label = Label(add_win, text='Reaction site :')
    mass_en = Entry(add_win, width=20)
    mass_label = Label(add_win, text='Cross linker mass :')
    short_en = Entry(add_win, width=20)
    short_label = Label(add_win, text='Short residue mass :')
    long_en = Entry(add_win, width=20)
    long_label = Label(add_win, text='Long residue mass :')

    linker_en.grid(row=0, column=1, padx=5, pady=5, sticky='nsew')
    linker_en.bind("<Enter>", on_enter_linker_name_en)
    linker_en.bind("<Leave>", on_cleave_add_linker_leave)
    linker_label.grid(row=0, column=0, padx=5, pady=5, sticky='nsew')
    site_en.grid(row=1, column=1, padx=5, pady=5, sticky='nsew')
    site_en.bind("<Enter>", on_enter_linker_site_en)
    site_en.bind("<Leave>", on_cleave_add_linker_leave)
    site_label.grid(row=1, column=0, padx=5, pady=5, sticky='nsew')
    mass_en.grid(row=2, column=1, padx=5, pady=5, sticky='nsew')
    mass_en.bind("<Enter>", on_enter_linker_mass_en)
    mass_en.bind("<Leave>", on_cleave_add_linker_leave)
    mass_label.grid(row=2, column=0, padx=5, pady=5, sticky='nsew')

    short_en.grid(row=3, column=1, padx=5, pady=5, sticky='nsew')
    short_en.bind("<Enter>", on_enter_linker_short_en)
    short_en.bind("<Leave>", on_cleave_add_linker_leave)
    short_label.grid(row=3, column=0, padx=5, pady=5, sticky='nsew')

    long_en.grid(row=4, column=1, padx=5, pady=5, sticky='nsew')
    long_en.bind("<Enter>", on_enter_linker_long_en)
    long_en.bind("<Leave>", on_cleave_add_linker_leave)
    long_label.grid(row=4, column=0, padx=5, pady=5, sticky='nsew')

    save_btn = Button(add_win, text='Save', command=save_cleave_add)
    save_btn.grid(row=5, column=0, columnspan=2, padx=5, pady=5, sticky='nsew')
    add_status = Label(add_win, text='Add cleavable cross linker', bd=1, relief=SUNKEN, anchor=E)
    add_status.grid(row=6, column=0, columnspan=2, ipady=4, sticky='sew')

    add_win.grid_rowconfigure(0, weight=1)
    add_win.grid_rowconfigure(1, weight=1)
    add_win.grid_rowconfigure(2, weight=1)
    add_win.grid_rowconfigure(3, weight=1)
    add_win.grid_rowconfigure(4, weight=1)
    add_win.grid_rowconfigure(5, weight=1)
    add_win.grid_rowconfigure(6, weight=1)
    add_win.grid_columnconfigure(0, weight=1)
    add_win.grid_columnconfigure(1, weight=1)


def edit_cleavelinker():
    """Cleavable module: edit linker"""
    global linker_en, site_en, mass_en, short_en, long_en, add_win, add_status
    add_win = Toplevel()
    add_win.title('Edit cleavable linker')
    add_win.iconbitmap('ECL.ico')
    add_win_w = 300
    add_win_h = 200
    add_win_width = add_win.winfo_screenwidth()
    add_win_height = add_win.winfo_screenheight()
    add_win_x = int(add_win_width / 2 - add_win_w / 2)
    add_win_y = int(add_win_height / 2 - add_win_h / 2)
    add_win.geometry('%dx%d+%d+%d' % (add_win_w, add_win_h, add_win_x, add_win_y))

    linker_en = Entry(add_win, width=20)
    linker_en.insert(0, linkerDict[cleavelinker_en.get()][0])
    linker_label = Label(add_win, text='Cross linker name :')
    site_en = Entry(add_win, width=20)
    site_en.insert(0, linkerDict[cleavelinker_en.get()][1])
    site_label = Label(add_win, text='Reaction site :')
    mass_en = Entry(add_win, width=20)
    mass_en.insert(0, linkerDict[cleavelinker_en.get()][2])
    mass_label = Label(add_win, text='Cross linker mass :')

    short_en = Entry(add_win, width=20)
    short_en.insert(0, linkerDict[cleavelinker_en.get()][3])
    short_en.bind("<Enter>", on_enter_linker_short_en)
    short_en.bind("<Leave>", on_cleave_edit_linker_leave)
    short_label = Label(add_win, text='Short residue mass :')

    long_en = Entry(add_win, width=20)
    long_en.insert(0, linkerDict[cleavelinker_en.get()][4])
    long_en.bind("<Enter>", on_enter_linker_long_en)
    long_en.bind("<Leave>", on_cleave_edit_linker_leave)
    long_label = Label(add_win, text='Long residue mass :')

    linker_en.grid(row=0, column=1, padx=5, pady=5, sticky='nsew')
    linker_en.bind("<Enter>", on_enter_linker_name_en)
    linker_en.bind("<Leave>", on_cleave_edit_linker_leave)
    linker_label.grid(row=0, column=0, padx=5, pady=5, sticky='nsew')
    site_en.grid(row=1, column=1, padx=5, pady=5, sticky='nsew')
    site_en.bind("<Enter>", on_enter_linker_site_en)
    site_en.bind("<Leave>", on_cleave_edit_linker_leave)
    site_label.grid(row=1, column=0, padx=5, pady=5, sticky='nsew')
    mass_en.grid(row=2, column=1, padx=5, pady=5, sticky='nsew')
    mass_en.bind("<Enter>", on_enter_linker_mass_en)
    mass_en.bind("<Leave>", on_cleave_edit_linker_leave)
    mass_label.grid(row=2, column=0, padx=5, pady=5, sticky='nsew')

    short_en.grid(row=3, column=1, padx=5, pady=5, sticky='nsew')
    short_label.grid(row=3, column=0, padx=5, pady=5, sticky='nsew')

    long_en.grid(row=4, column=1, padx=5, pady=5, sticky='nsew')
    long_label.grid(row=4, column=0, padx=5, pady=5, sticky='nsew')

    save_btn = Button(add_win, text='Save', command=save_cleave_edit)
    save_btn.grid(row=5, column=0, columnspan=2, padx=5, pady=5, sticky='nsew')
    add_status = Label(add_win, text='Edit cleavable cross linker', bd=1, relief=SUNKEN, anchor=E)
    add_status.grid(row=6, column=0, columnspan=2, ipady=4, sticky='sew')

    add_win.grid_rowconfigure(0, weight=1)
    add_win.grid_rowconfigure(1, weight=1)
    add_win.grid_rowconfigure(2, weight=1)
    add_win.grid_rowconfigure(3, weight=1)
    add_win.grid_rowconfigure(4, weight=1)
    add_win.grid_rowconfigure(5, weight=1)
    add_win.grid_rowconfigure(6, weight=1)
    add_win.grid_columnconfigure(0, weight=1)
    add_win.grid_columnconfigure(1, weight=1)


def refresh_cleavelinker():
    """Cleavable module: refresh linker database"""
    linkerDict.clear()
    conn = sqlite3.connect('paramsDB.db')
    c = conn.cursor()
    c.execute("SELECT *, oid FROM cleavelinker")
    for ele in c.fetchall():
        linkerDict[ele[0]] = ele
    conn.commit()
    conn.close()
    cleavelinker_en['values'] = list(linkerDict.keys())
    cleavelinker_en.current(0)


def delete_cleavelinker():
    """Cleavable module: delete linker"""
    response = messagebox.askyesno("Delete cross linker",
                                   "Delete linker {} in the database?".format(cleavelinker_en.get()))
    if response == 1:
        conn = sqlite3.connect('paramsDB.db')
        c = conn.cursor()
        c.execute("DELETE from cleavelinker WHERE oid= " + str(linkerDict[cleavelinker_en.get()][-1]))

        conn.commit()
        conn.close()


def save_mod_add():
    """Noncleavable module: save added modification database"""
    if not mod_name_en.get().islower():
        messagebox.showwarning('Invalid input', 'Invalid mods. name')
        return True
    if len(mod_site_en.get()) != 1 or not mod_site_en.get().isupper():
        messagebox.showwarning('Invalid input', 'Invalid reaction site')
        return True
    try:
        float(mod_mass_en.get())
    except ValueError:
        messagebox.showwarning('Invalid input', 'Invalid mods. mass')
        return True

    conn = sqlite3.connect('paramsDB.db')
    c = conn.cursor()
    c.execute("INSERT INTO nonlinkermod VALUES (:name, :site, :mass)",
              {
                  'name': mod_name_en.get(),
                  'site': mod_site_en.get(),
                  'mass': mod_mass_en.get()
              })
    conn.commit()
    conn.close()
    add_win.destroy()


def add_mod_nonlinker():
    """Noncleavable module: add modification database"""
    global mod_name_en, mod_site_en, mod_mass_en, add_win, add_status
    add_win = Toplevel()
    add_win.title('Add modification')
    add_win.iconbitmap('ECL.ico')
    add_win_w = 300
    add_win_h = 150
    add_win_width = add_win.winfo_screenwidth()
    add_win_height = add_win.winfo_screenheight()
    add_win_x = int(add_win_width / 2 - add_win_w / 2)
    add_win_y = int(add_win_height / 2 - add_win_h / 2)
    add_win.geometry('%dx%d+%d+%d' % (add_win_w, add_win_h, add_win_x, add_win_y))

    mod_name_en = Entry(add_win, width=20)
    mod_name_label = Label(add_win, text='Modification name :')
    mod_site_en = Entry(add_win, width=20)
    mod_site_label = Label(add_win, text='Reaction site :')
    mod_mass_en = Entry(add_win, width=20)
    mod_mass_label = Label(add_win, text='Modification mass :')
    mod_name_en.grid(row=0, column=1, padx=5, pady=5, sticky='nsew')
    mod_name_en.bind("<Enter>", on_enter_mod_name_en)
    mod_name_en.bind("<Leave>", on_non_add_mod_leave)
    mod_name_label.grid(row=0, column=0, padx=5, pady=5, sticky='nsew')
    mod_site_en.grid(row=1, column=1, padx=5, pady=5, sticky='nsew')
    mod_site_en.bind("<Enter>", on_enter_mod_site_en)
    mod_site_en.bind("<Leave>", on_non_add_mod_leave)
    mod_site_label.grid(row=1, column=0, padx=5, pady=5, sticky='nsew')
    mod_mass_en.grid(row=2, column=1, padx=5, pady=5, sticky='nsew')
    mod_mass_en.bind("<Enter>", on_enter_mod_mass_en)
    mod_mass_en.bind("<Leave>", on_non_add_mod_leave)
    mod_mass_label.grid(row=2, column=0, padx=5, pady=5, sticky='nsew')
    save_btn = Button(add_win, text='Save', command=save_mod_add)
    save_btn.grid(row=3, column=0, columnspan=2, padx=5, pady=5, sticky='nsew')
    add_status = Label(add_win, text='Add modification', bd=1, relief=SUNKEN, anchor=E)
    add_status.grid(row=4, column=0, columnspan=2, ipady=4, sticky='sew')

    add_win.grid_rowconfigure(0, weight=1)
    add_win.grid_rowconfigure(1, weight=1)
    add_win.grid_rowconfigure(2, weight=1)
    add_win.grid_rowconfigure(3, weight=1)
    add_win.grid_rowconfigure(4, weight=1)
    add_win.grid_columnconfigure(0, weight=1)
    add_win.grid_columnconfigure(1, weight=1)


def refresh_mod_nonlinker():
    """Noncleavable module: update modification database"""
    modDict.clear()
    conn = sqlite3.connect('paramsDB.db')
    c = conn.cursor()
    c.execute("SELECT *, oid FROM nonlinkermod")
    # print(c.fetchall())
    for ele in c.fetchall():
        modDict[ele[0] + '[' + ele[1] + ']'] = ele
    conn.commit()
    conn.close()
    # print(modDict.items())
    mod_en.delete(0, END)
    for ele in modDict.keys():
        mod_en.insert(END, ele)


def delete_mod_nonlinker():
    """Noncleavable module: delete modification database"""
    if mod_en.curselection():
        response = messagebox.askyesno("Delete modification",
                                       "Delete mods. {} in the database?".format(mod_en.get(ANCHOR)))
        if response == 1:
            conn = sqlite3.connect('paramsDB.db')
            c = conn.cursor()
            c.execute("DELETE from nonlinkermod WHERE oid= " + str(modDict[mod_en.get(ANCHOR)][-1]))

            conn.commit()
            conn.close()


def save_mod_edit():
    """Noncleavable module: save edited modification database"""
    if not mod_name_en.get().islower():
        messagebox.showwarning('Invalid input', 'Invalid mods. name')
        return True
    if len(mod_site_en.get()) != 1 or not mod_site_en.get().isupper():
        messagebox.showwarning('Invalid input', 'Invalid reaction site')
        return True
    try:
        float(mod_mass_en.get())
    except ValueError:
        messagebox.showwarning('Invalid input', 'Invalid mods. mass')
        return True

    conn = sqlite3.connect('paramsDB.db')
    c = conn.cursor()
    c.execute("""UPDATE nonlinkermod SET
        name = :name,
        site = :site,
        mass = :mass
        WHERE oid = :oid""",
              {'name': mod_name_en.get(),
               'site': mod_site_en.get(),
               'mass': mod_mass_en.get(),
               'oid': modDict[mod_en.get(ANCHOR)][-1]})
    conn.commit()
    conn.close()
    add_win.destroy()


def edit_mod_nonlinker():
    """Noncleavable module: edit modification database"""
    global mod_name_en, mod_site_en, mod_mass_en, add_win, add_status
    if mod_en.curselection():
        add_win = Toplevel()
        add_win.title('Edit modification')
        add_win.iconbitmap('ECL.ico')
        add_win_w = 300
        add_win_h = 150
        add_win_width = add_win.winfo_screenwidth()
        add_win_height = add_win.winfo_screenheight()
        add_win_x = int(add_win_width / 2 - add_win_w / 2)
        add_win_y = int(add_win_height / 2 - add_win_h / 2)
        add_win.geometry('%dx%d+%d+%d' % (add_win_w, add_win_h, add_win_x, add_win_y))

        mod_name_en = Entry(add_win, width=20)
        mod_name_en.insert(0, modDict[mod_en.get(ANCHOR)][0])
        mod_name_label = Label(add_win, text='Modification name :')
        mod_site_en = Entry(add_win, width=20)
        mod_site_en.insert(0, modDict[mod_en.get(ANCHOR)][1])
        mod_site_label = Label(add_win, text='Reaction site :')
        mod_mass_en = Entry(add_win, width=20)
        mod_mass_en.insert(0, modDict[mod_en.get(ANCHOR)][2])
        mod_mass_label = Label(add_win, text='Modification mass :')
        mod_name_en.grid(row=0, column=1, padx=5, pady=5, sticky='nsew')
        mod_name_en.bind("<Enter>", on_enter_mod_name_en)
        mod_name_en.bind("<Leave>", on_non_edit_mod_leave)
        mod_name_label.grid(row=0, column=0, padx=5, pady=5, sticky='nsew')
        mod_site_en.grid(row=1, column=1, padx=5, pady=5, sticky='nsew')
        mod_site_en.bind("<Enter>", on_enter_mod_site_en)
        mod_site_en.bind("<Leave>", on_non_edit_mod_leave)
        mod_site_label.grid(row=1, column=0, padx=5, pady=5, sticky='nsew')
        mod_mass_en.grid(row=2, column=1, padx=5, pady=5, sticky='nsew')
        mod_mass_en.bind("<Enter>", on_enter_mod_mass_en)
        mod_mass_en.bind("<Leave>", on_non_edit_mod_leave)
        mod_mass_label.grid(row=2, column=0, padx=5, pady=5, sticky='nsew')
        save_btn = Button(add_win, text='Save', command=save_mod_edit)
        save_btn.grid(row=3, column=0, columnspan=2, padx=5, pady=5, sticky='nsew')
        add_status = Label(add_win, text='Edit modification', bd=1, relief=SUNKEN, anchor=E)
        add_status.grid(row=4, column=0, columnspan=2, ipady=4, sticky='sew')

        add_win.grid_rowconfigure(0, weight=1)
        add_win.grid_rowconfigure(1, weight=1)
        add_win.grid_rowconfigure(2, weight=1)
        add_win.grid_rowconfigure(3, weight=1)
        add_win.grid_rowconfigure(4, weight=1)
        add_win.grid_columnconfigure(0, weight=1)
        add_win.grid_columnconfigure(1, weight=1)


def save_cleave_mod_add():
    """Cleavable module: save added modification database"""
    if not mod_name_en.get().islower():
        messagebox.showwarning('Invalid input', 'Invalid mods. name')
        return True
    if len(mod_site_en.get()) != 1 or not mod_site_en.get().isupper():
        messagebox.showwarning('Invalid input', 'Invalid reaction site')
        return True
    try:
        float(mod_mass_en.get())
    except ValueError:
        messagebox.showwarning('Invalid input', 'Invalid mods. mass')
        return True

    conn = sqlite3.connect('paramsDB.db')
    c = conn.cursor()
    c.execute("INSERT INTO cleavelinkermod VALUES (:name, :site, :mass)",
              {
                  'name': mod_name_en.get(),
                  'site': mod_site_en.get(),
                  'mass': mod_mass_en.get()
              })
    conn.commit()
    conn.close()
    add_win.destroy()


def add_mod_cleavelinker():
    """Cleavable module: add modification database"""
    global mod_name_en, mod_site_en, mod_mass_en, add_win, add_status
    add_win = Toplevel()
    add_win.title('Add modification')
    add_win.iconbitmap('ECL.ico')
    add_win_w = 300
    add_win_h = 150
    add_win_width = add_win.winfo_screenwidth()
    add_win_height = add_win.winfo_screenheight()
    add_win_x = int(add_win_width / 2 - add_win_w / 2)
    add_win_y = int(add_win_height / 2 - add_win_h / 2)
    add_win.geometry('%dx%d+%d+%d' % (add_win_w, add_win_h, add_win_x, add_win_y))

    mod_name_en = Entry(add_win, width=20)
    mod_name_label = Label(add_win, text='Modification name :')
    mod_site_en = Entry(add_win, width=20)
    mod_site_label = Label(add_win, text='Reaction site :')
    mod_mass_en = Entry(add_win, width=20)
    mod_mass_label = Label(add_win, text='Modification mass :')
    mod_name_en.grid(row=0, column=1, padx=5, pady=5, sticky='nsew')
    mod_name_en.bind("<Enter>", on_enter_mod_name_en)
    mod_name_en.bind("<Leave>", on_non_add_mod_leave)
    mod_name_label.grid(row=0, column=0, padx=5, pady=5, sticky='nsew')
    mod_site_en.grid(row=1, column=1, padx=5, pady=5, sticky='nsew')
    mod_site_en.bind("<Enter>", on_enter_mod_site_en)
    mod_site_en.bind("<Leave>", on_non_add_mod_leave)
    mod_site_label.grid(row=1, column=0, padx=5, pady=5, sticky='nsew')
    mod_mass_en.grid(row=2, column=1, padx=5, pady=5, sticky='nsew')
    mod_mass_en.bind("<Enter>", on_enter_mod_mass_en)
    mod_mass_en.bind("<Leave>", on_non_add_mod_leave)
    mod_mass_label.grid(row=2, column=0, padx=5, pady=5, sticky='nsew')
    save_btn = Button(add_win, text='Save', command=save_cleave_mod_add)
    save_btn.grid(row=3, column=0, columnspan=2, padx=5, pady=5, sticky='nsew')
    add_status = Label(add_win, text='Add modification', bd=1, relief=SUNKEN, anchor=E)
    add_status.grid(row=4, column=0, columnspan=2, ipady=4, sticky='sew')

    add_win.grid_rowconfigure(0, weight=1)
    add_win.grid_rowconfigure(1, weight=1)
    add_win.grid_rowconfigure(2, weight=1)
    add_win.grid_rowconfigure(3, weight=1)
    add_win.grid_rowconfigure(4, weight=1)
    add_win.grid_columnconfigure(0, weight=1)
    add_win.grid_columnconfigure(1, weight=1)


def refresh_mod_cleavelinker():
    """Cleavable module: update modification database"""
    modDict.clear()
    conn = sqlite3.connect('paramsDB.db')
    c = conn.cursor()
    c.execute("SELECT *, oid FROM cleavelinkermod")
    for ele in c.fetchall():
        modDict[ele[0] + '[' + ele[1] + ']'] = ele
    conn.commit()
    conn.close()
    mod_en.delete(0, END)
    for ele in modDict.keys():
        mod_en.insert(END, ele)


def delete_mod_cleavelinker():
    """Cleavable module: delete modification database"""
    if mod_en.curselection():
        response = messagebox.askyesno("Delete modification",
                                       "Delete mods. {} in the database?".format(mod_en.get(ANCHOR)))
        if response == 1:
            conn = sqlite3.connect('paramsDB.db')
            c = conn.cursor()
            c.execute("DELETE from cleavelinkermod WHERE oid= " + str(modDict[mod_en.get(ANCHOR)][-1]))

            conn.commit()
            conn.close()


def save_cleave_mod_edit():
    """Cleavable module: save edited modification database"""
    if not mod_name_en.get().islower():
        messagebox.showwarning('Invalid input', 'Invalid mods. name')
        return True
    if len(mod_site_en.get()) != 1 or not mod_site_en.get().isupper():
        messagebox.showwarning('Invalid input', 'Invalid reaction site')
        return True
    try:
        float(mod_mass_en.get())
    except ValueError:
        messagebox.showwarning('Invalid input', 'Invalid mods. mass')
        return True

    conn = sqlite3.connect('paramsDB.db')
    c = conn.cursor()
    c.execute("""UPDATE cleavelinkermod SET
        name = :name,
        site = :site,
        mass = :mass
        WHERE oid = :oid""",
              {'name': mod_name_en.get(),
               'site': mod_site_en.get(),
               'mass': mod_mass_en.get(),
               'oid': modDict[mod_en.get(ANCHOR)][-1]})
    conn.commit()
    conn.close()
    add_win.destroy()


def edit_mod_cleavelinker():
    """Cleavable module: edit modification database"""
    global mod_name_en, mod_site_en, mod_mass_en, add_win, add_status
    if mod_en.curselection():
        add_win = Toplevel()
        add_win.title('Edit modification')
        add_win.iconbitmap('ECL.ico')
        add_win_w = 300
        add_win_h = 150
        add_win_width = add_win.winfo_screenwidth()
        add_win_height = add_win.winfo_screenheight()
        add_win_x = int(add_win_width / 2 - add_win_w / 2)
        add_win_y = int(add_win_height / 2 - add_win_h / 2)
        add_win.geometry('%dx%d+%d+%d' % (add_win_w, add_win_h, add_win_x, add_win_y))

        mod_name_en = Entry(add_win, width=20)
        mod_name_en.insert(0, modDict[mod_en.get(ANCHOR)][0])
        mod_name_label = Label(add_win, text='Modification name :')
        mod_site_en = Entry(add_win, width=20)
        mod_site_en.insert(0, modDict[mod_en.get(ANCHOR)][1])
        mod_site_label = Label(add_win, text='Reaction site :')
        mod_mass_en = Entry(add_win, width=20)
        mod_mass_en.insert(0, modDict[mod_en.get(ANCHOR)][2])
        mod_mass_label = Label(add_win, text='Modification mass :')
        mod_name_en.grid(row=0, column=1, padx=5, pady=5, sticky='nsew')
        mod_name_en.bind("<Enter>", on_enter_mod_name_en)
        mod_name_en.bind("<Leave>", on_non_edit_mod_leave)
        mod_name_label.grid(row=0, column=0, padx=5, pady=5, sticky='nsew')
        mod_site_en.grid(row=1, column=1, padx=5, pady=5, sticky='nsew')
        mod_site_en.bind("<Enter>", on_enter_mod_site_en)
        mod_site_en.bind("<Leave>", on_non_edit_mod_leave)
        mod_site_label.grid(row=1, column=0, padx=5, pady=5, sticky='nsew')
        mod_mass_en.grid(row=2, column=1, padx=5, pady=5, sticky='nsew')
        mod_mass_en.bind("<Enter>", on_enter_mod_mass_en)
        mod_mass_en.bind("<Leave>", on_non_edit_mod_leave)
        mod_mass_label.grid(row=2, column=0, padx=5, pady=5, sticky='nsew')
        save_btn = Button(add_win, text='Save', command=save_cleave_mod_edit)
        save_btn.grid(row=3, column=0, columnspan=2, padx=5, pady=5, sticky='nsew')
        add_status = Label(add_win, text='Edit modification', bd=1, relief=SUNKEN, anchor=E)
        add_status.grid(row=4, column=0, columnspan=2, ipady=4, sticky='sew')

        add_win.grid_rowconfigure(0, weight=1)
        add_win.grid_rowconfigure(1, weight=1)
        add_win.grid_rowconfigure(2, weight=1)
        add_win.grid_rowconfigure(3, weight=1)
        add_win.grid_rowconfigure(4, weight=1)
        add_win.grid_columnconfigure(0, weight=1)
        add_win.grid_columnconfigure(1, weight=1)


def on_cleave_leave(e):
    """Cleaevable module: when leave button"""
    e.widget['background'] = 'SystemButtonFace'
    non_status.config(text='Cleavable cross-linker searching module')


def on_cleave_add_linker_leave(e):
    """Cleavable module: leave linker frame"""
    add_status.config(text='Add cleavable cross linker')


def on_enter_linker_short_en(e):
    """Cleavable module: enter short linker mass"""
    add_status.config(text='Input short mono (residual) mass')


def on_enter_linker_long_en(e):
    """Cleavable module: enter long linker mass"""
    add_status.config(text='Input long mono (residual) mass')


def on_cleave_edit_linker_leave(e):
    """Cleavable module: leave edit linker mass"""
    add_status.config(text='Edit cleavable cross linker')


'''################################################ Noncleavable frame ##############################################'''


def save_add():
    """Noncleavable module: save added linker"""
    if not linker_en.get():
        messagebox.showwarning('Invalid input', 'No linker name')
        return True
    if len(site_en.get()) != 1 or not site_en.get().isupper():
        messagebox.showwarning('Invalid input', 'Invalid link site')
        return True
    try:
        float(mass_en.get())
    except ValueError:
        messagebox.showwarning('Invalid input', 'Invalid linker mass')
        return True

    conn = sqlite3.connect('paramsDB.db')
    c = conn.cursor()
    c.execute("INSERT INTO nonlinker VALUES (:linker, :site, :mass)",
              {
                  'linker': linker_en.get(),
                  'site': site_en.get(),
                  'mass': mass_en.get()
              })
    conn.commit()
    conn.close()
    add_win.destroy()


def save_edit():
    """Noncleavable module: save edited linker"""
    if not linker_en.get():
        messagebox.showwarning('Invalid input', 'No linker name')
        return True
    if len(site_en.get()) != 1 or not site_en.get().isupper():
        messagebox.showwarning('Invalid input', 'Invalid link site')
        return True
    try:
        float(mass_en.get())
    except ValueError:
        messagebox.showwarning('Invalid input', 'Invalid linker mass')
        return True

    conn = sqlite3.connect('paramsDB.db')
    c = conn.cursor()
    c.execute("""UPDATE nonlinker SET
        linker = :linker,
        site = :site,
        mass = :mass
        WHERE oid = :oid""",
              {'linker': linker_en.get(),
               'site': site_en.get(),
               'mass': mass_en.get(),
               'oid': linkerDict[nonlinker_en.get()][3]})
    conn.commit()
    conn.close()
    add_win.destroy()


def add_nonlinker():
    """Noncleavable module: add linker"""
    global linker_en, site_en, mass_en, add_win, add_status
    add_win = Toplevel()
    add_win.title('Add non-cleavable linker')
    add_win.iconbitmap('ECL.ico')
    add_win_w = 300
    add_win_h = 150
    add_win_width = add_win.winfo_screenwidth()
    add_win_height = add_win.winfo_screenheight()
    add_win_x = int(add_win_width / 2 - add_win_w / 2)
    add_win_y = int(add_win_height / 2 - add_win_h / 2)
    add_win.geometry('%dx%d+%d+%d' % (add_win_w, add_win_h, add_win_x, add_win_y))

    linker_en = Entry(add_win, width=20)
    linker_label = Label(add_win, text='Cross linker name :')
    site_en = Entry(add_win, width=20)
    site_label = Label(add_win, text='Reaction site :')
    mass_en = Entry(add_win, width=20)
    mass_label = Label(add_win, text='Cross linker mass :')
    linker_en.grid(row=0, column=1, padx=5, pady=5, sticky='nsew')
    linker_en.bind("<Enter>", on_enter_linker_name_en)
    linker_en.bind("<Leave>", on_non_add_linker_leave)
    linker_label.grid(row=0, column=0, padx=5, pady=5, sticky='nsew')
    site_en.grid(row=1, column=1, padx=5, pady=5, sticky='nsew')
    site_en.bind("<Enter>", on_enter_linker_site_en)
    site_en.bind("<Leave>", on_non_add_linker_leave)
    site_label.grid(row=1, column=0, padx=5, pady=5, sticky='nsew')
    mass_en.grid(row=2, column=1, padx=5, pady=5, sticky='nsew')
    mass_en.bind("<Enter>", on_enter_linker_mass_en)
    mass_en.bind("<Leave>", on_non_add_linker_leave)
    mass_label.grid(row=2, column=0, padx=5, pady=5, sticky='nsew')
    save_btn = Button(add_win, text='Save', command=save_add)
    save_btn.grid(row=3, column=0, columnspan=2, padx=5, pady=5, sticky='nsew')
    add_status = Label(add_win, text='Add non-cleavable cross linker', bd=1, relief=SUNKEN, anchor=E)
    add_status.grid(row=4, column=0, columnspan=2, ipady=4, sticky='sew')

    add_win.grid_rowconfigure(0, weight=1)
    add_win.grid_rowconfigure(1, weight=1)
    add_win.grid_rowconfigure(2, weight=1)
    add_win.grid_rowconfigure(3, weight=1)
    add_win.grid_rowconfigure(4, weight=1)
    add_win.grid_columnconfigure(0, weight=1)
    add_win.grid_columnconfigure(1, weight=1)


def edit_nonlinker():
    """Noncleavable module: edit linker"""
    global linker_en, site_en, mass_en, add_win, add_status
    add_win = Toplevel()
    add_win.title('Edit non-cleavable linker')
    add_win.iconbitmap('ECL.ico')
    add_win_w = 300
    add_win_h = 150
    add_win_width = add_win.winfo_screenwidth()
    add_win_height = add_win.winfo_screenheight()
    add_win_x = int(add_win_width / 2 - add_win_w / 2)
    add_win_y = int(add_win_height / 2 - add_win_h / 2)
    add_win.geometry('%dx%d+%d+%d' % (add_win_w, add_win_h, add_win_x, add_win_y))

    linker_en = Entry(add_win, width=20)
    linker_en.insert(0, linkerDict[nonlinker_en.get()][0])
    linker_label = Label(add_win, text='Cross linker name :')
    site_en = Entry(add_win, width=20)
    site_en.insert(0, linkerDict[nonlinker_en.get()][1])
    site_label = Label(add_win, text='Reaction site :')
    mass_en = Entry(add_win, width=20)
    mass_en.insert(0, linkerDict[nonlinker_en.get()][2])
    mass_label = Label(add_win, text='Cross linker mass :')
    linker_en.grid(row=0, column=1, padx=5, pady=5, sticky='nsew')
    linker_en.bind("<Enter>", on_enter_linker_name_en)
    linker_en.bind("<Leave>", on_non_edit_linker_leave)
    linker_label.grid(row=0, column=0, padx=5, pady=5, sticky='nsew')
    site_en.grid(row=1, column=1, padx=5, pady=5, sticky='nsew')
    site_en.bind("<Enter>", on_enter_linker_site_en)
    site_en.bind("<Leave>", on_non_edit_linker_leave)
    site_label.grid(row=1, column=0, padx=5, pady=5, sticky='nsew')
    mass_en.grid(row=2, column=1, padx=5, pady=5, sticky='nsew')
    mass_en.bind("<Enter>", on_enter_linker_mass_en)
    mass_en.bind("<Leave>", on_non_edit_linker_leave)
    mass_label.grid(row=2, column=0, padx=5, pady=5, sticky='nsew')
    save_btn = Button(add_win, text='Save', command=save_edit)
    save_btn.grid(row=3, column=0, columnspan=2, padx=5, pady=5, sticky='nsew')
    add_status = Label(add_win, text='Edit non-cleavable cross linker', bd=1, relief=SUNKEN, anchor=E)
    add_status.grid(row=4, column=0, columnspan=2, ipady=4, sticky='sew')

    add_win.grid_rowconfigure(0, weight=1)
    add_win.grid_rowconfigure(1, weight=1)
    add_win.grid_rowconfigure(2, weight=1)
    add_win.grid_rowconfigure(3, weight=1)
    add_win.grid_rowconfigure(4, weight=1)
    add_win.grid_columnconfigure(0, weight=1)
    add_win.grid_columnconfigure(1, weight=1)


def refresh_nonlinker():
    """Noncleavable module: refresh linker database"""
    linkerDict.clear()
    conn = sqlite3.connect('paramsDB.db')
    c = conn.cursor()
    c.execute("SELECT *, oid FROM nonlinker")
    for ele in c.fetchall():
        linkerDict[ele[0]] = ele
    conn.commit()
    conn.close()
    nonlinker_en['values'] = list(linkerDict.keys())
    nonlinker_en.current(0)


def delete_nonlinker():
    """Noncleavable module: delete linker"""
    response = messagebox.askyesno("Delete cross linker",
                                   "Delete linker {} in the database?".format(nonlinker_en.get()))
    if response == 1:
        conn = sqlite3.connect('paramsDB.db')
        c = conn.cursor()
        c.execute("DELETE from nonlinker WHERE oid= " + str(linkerDict[nonlinker_en.get()][-1]))

        conn.commit()
        conn.close()


def non_params_print():
    """Noncleavable module: print parameters"""
    print_en.configure(state=NORMAL)
    print_en.delete(1.0, END)
    print_en.insert(END, "{}Parameters{}\n\n".format('*'*50, '*'*50))
    print_en.insert(END, "Data :              ")
    print_en.insert(END, infile[0] + '\n')
    for ele in infile[1:]:
        print_en.insert(END,'                    ' + ele + '\n')
    print_en.insert(END, 'FASTA file :        ' + db_en.get() + '\n')
    print_en.insert(END, "Output:             ")
    print_en.insert(END, infile[0].replace('.mzXML', '_final.csv') + '\n')
    for ele in infile[1:]:
        print_en.insert(END, '                    ' + ele.replace('.mzXML', '_final.csv') + '\n')
    print_en.insert(END, 'Threads :           ' + core_en.get() + '\n')
    print_en.insert(END, 'Peptide mass :      ' + '[' + lower_en.get() + ', ' + higher_en.get() + ']' + '\n')
    print_en.insert(END, 'Digestion :         ' + parse_en.get() + '\n')
    print_en.insert(END, 'Miss cleavage :     ' + miss_en.get() + '\n')
    print_en.insert(END, 'Ms1 tolerance :     ' + ms1tol_en.get() + ' ppm\n')
    print_en.insert(END, 'Ms2 tolerance :     ' + ms2tol_en.get() + ' Da\n')
    print_en.insert(END, 'Cross linker :      ' + nonlinker_en.get() + '\n')
    print_en.insert(END, 'Link site :         ' + linkerDict[nonlinker_en.get()][1] + '\n')
    print_en.insert(END, 'Linker mass :       ' + str(linkerDict[nonlinker_en.get()][2]) + '\n')
    if fixed_en.get(0, END):
        print_en.insert(END, 'Fixed mods. :       ' + fixed_en.get(0, END)[0] + ' ' + str(modDict[fixed_en.get(0, END)[0]][2]) + '\n')
        for ele in fixed_en.get(0, END)[1:]:
            print_en.insert(END, '                    ' + ele + ' ' + str(modDict[ele][2]) + '\n')
    if var_en.get(0, END):
        print_en.insert(END, 'Var mods. :         ' + var_en.get(0, END)[0] + ' ' + str(modDict[var_en.get(0, END)[0]][2]) + '\n')
        for ele in var_en.get(0, END)[1:]:
            print_en.insert(END, '                    ' + ele + ' ' + str(modDict[ele][2]) + '\n')
    print_en.insert(END, 'FDR threshold :     ' + fdr_en.get() + ' %\n\n')
    print_en.insert(END, "{}Program running{}\n\n".format('*' *47, '*' * 48))


def sub_xolik():
    """Noncleavable: run data set"""
    global p1, p2, flag1
    if fixed_en.get(0, END):
        fix_in = modDict[fixed_en.get(0, END)[0]][1] + '+' + str(modDict[fixed_en.get(0, END)[0]][2]) + '&' + \
                 modDict[fixed_en.get(0, END)[0]][0]
        for ele in fixed_en.get(0, END)[1:]:
            fix_in += ' : ' + modDict[ele][1] + '+' + str(modDict[ele][2]) + '&' + modDict[ele][0]
    else:
        fix_in = ''
    if var_en.get(0, END):
        var_in = modDict[var_en.get(0, END)[0]][1] + '+' + str(modDict[var_en.get(0, END)[0]][2]) + '&' + \
                 modDict[var_en.get(0, END)[0]][0]
        for ele in var_en.get(0, END)[1:]:
            var_in += ' : ' + modDict[ele][1] + '+' + str(modDict[ele][2]) + '&' + modDict[ele][0]
    else:
        var_in = ''
    flag1 = False
    for file_input in infile:
        print_en.insert(END, '>> Analyzing file {}\n'.format(file_input))
        outputname = file_input.replace('.mzXML', '.csv')
        command = 'Xolik.exe -d "{}" -s "{}" -o "{}" --miss {} ' \
                  '--min {} --max {} --xlsite {} --xlmass {}' \
                  ' --ms1tol {} --ms2tol {} --parallel' \
                  ' --thread {} --varmod "{}" --fixmod "{}"'.format(db_en.get(), file_input, outputname,
                                                                    miss_en.get(), lower_en.get(), higher_en.get(),
                                                                    linkerDict[nonlinker_en.get()][1],
                                                                    str(linkerDict[nonlinker_en.get()][2]),
                                                                    ms1tol_en.get(), ms2tol_en.get(),
                                                                    str(core_en.get()), var_in, fix_in)
        p1 = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              stdin=subprocess.PIPE, shell=False)

        while p1.poll() is None:
            for line in iter(p1.stdout.readline, b""):
                print_en.insert(END, line.decode("utf-8"))
                print_en.see(END)

        p1.wait()
        if flag1:
            return True
        print_en.insert(END, '>> Protein feedback processing...\n')
        command_filter = 'python noncleavable_splitctrl.py "{}" "{}" {}'.format(outputname, db_en.get(),
                                                                                float(fdr_en.get())*0.01)

        p2 = subprocess.call(command_filter, shell=False)

    print_en.insert(END, '>> End at {}\n\n'.format(time.ctime()))
    print_en.insert(END, "{}Program finished{}\n".format('*' * 47, '*' * 47))

    print_en.configure(state=DISABLED)
    stop_btn['state'] = DISABLED
    messagebox.showinfo('Task', 'Task finished!')


def run_xolik():
    """Noncleavable module: check input params and run program"""
    if not infile:
        messagebox.showwarning('Invalid parameter settings', 'No input data file(s)')
        return True

    if not db_en.get():
        messagebox.showwarning('Invalid parameter settings', 'No input FASTA file')
        return True

    try:
        if float(lower_en.get()) >= float(higher_en.get()):
            messagebox.showwarning('Invalid parameter settings', 'Peptide lower mass is larger than higher mass')
            return True
    except ValueError:
        messagebox.showwarning('Invalid parameter settings', 'Invalid peptide mass')
        return True

    try:
        if float(ms1tol_en.get()) <=0 or float(ms2tol_en.get()) <= 0:
            messagebox.showwarning('Invalid parameter settings', 'Mass tolerance is smaller than 0')
            return True
    except ValueError:
        messagebox.showwarning('Invalid parameter settings', 'Invalid mass tolerance')
        return True

    try:
        if float(fdr_en.get()) <= 0:
            messagebox.showwarning('Invalid parameter settings', 'FDR is less than 0')
            return True
    except ValueError:
        messagebox.showwarning('Invalid parameter settings', 'Invalid FDR setting')
        return True

    stop_btn['state'] = NORMAL
    print_en.configure(state=NORMAL)
    non_params_print()

    threading.Thread(target=sub_xolik).start()


def noncleav():
    """Noncleavable module: frame"""
    global non_win
    global img_back, img_start, img_stop
    global data_en, db_en, core_en, infile
    global lower_en, higher_en, ms1tol_en, ms2tol_en, parse_en, miss_en, fdr_en
    global nonlinker_en, linkerDict
    global fixed_en, var_en, mod_en, modDict
    global print_en, stop_btn
    global non_status
    noncleavable['state'] = DISABLED
    cleavable['state'] = DISABLED
    root.withdraw()

    non_win = Toplevel()
    non_win.title('Non-cleavable searching module')
    non_win.iconbitmap('ECL.ico')
    non_win_w = 900
    non_win_h = 700
    non_win_width = non_win.winfo_screenwidth()
    non_win_height = non_win.winfo_screenheight()
    non_win_x = int(non_win_width / 2 - non_win_w / 2)
    non_win_y = int(non_win_height / 2 - non_win_h / 2)
    non_win.geometry('%dx%d+%d+%d' % (non_win_w, non_win_h, non_win_x, non_win_y))
    non_win.protocol("WM_DELETE_WINDOW", back2normal1)
    '''Input frame'''
    input_frame = LabelFrame(non_win, text='Input')
    input_frame.grid(row=0, column=0, rowspan=3, columnspan=2, padx=5, pady=(10,5), sticky='nsew')
    input_frame.grid_propagate(False)

    data_label = Label(input_frame, text="Data file(s)")
    data_label.grid(row=0, column=0, rowspan=4, sticky="nw")

    data_en_scroll = Scrollbar(input_frame, orient='horizontal')
    data_en_scroll.grid(row=3, column=1, sticky="nsew")
    data_en = Text(input_frame, state=DISABLED, width=10, height=5, wrap=NONE, font=('Arial', 8), xscrollcommand=data_en_scroll.set)
    data_en.grid(row=0, column=1, rowspan=3, sticky="nsew")
    data_en_scroll.config(command=data_en.xview)

    infile = tuple()
    data_btn = Button(input_frame, text="Browse", command=browse)
    data_btn.grid(row=0, column=2, rowspan=4, padx=5, sticky="new")
    data_btn.bind("<Enter>", on_enter_data_btn)
    data_btn.bind("<Leave>", on_non_leave)

    db_label = Label(input_frame, text="Database")
    db_label.grid(row=4, column=0, sticky="nsew")
    db_en = Entry(input_frame, font=('Arial', 8), state=DISABLED)
    db_en.grid(row=4, column=1, sticky="nsew")
    db_btn = Button(input_frame, text="Browse", command=browsefasta)
    db_btn.grid(row=4, column=2, padx=5, sticky="nsew")
    db_btn.bind("<Enter>", on_enter_db_btn)
    db_btn.bind("<Leave>", on_non_leave)



    core_label = Label(input_frame, text="# threads")
    core_label.grid(row=5, column=0, sticky="nsew")
    core_en = ttk.Combobox(input_frame, width=4, state="readonly")
    core_en['values'] = (1, 2, 4, 8, 16)
    core_en.current(2)
    core_en.grid(row=5, column=1, pady=5, sticky = "nsw")

    input_frame.grid_columnconfigure(1, weight=6)
    input_frame.grid_columnconfigure(2, weight=1)

    input_frame.grid_rowconfigure(0, weight=1)
    input_frame.grid_rowconfigure(1, weight=1)
    input_frame.grid_rowconfigure(2, weight=1)

    '''searching frame'''
    search_frame = LabelFrame(non_win, text='Parameters')
    search_frame.grid(row=3, column=0, rowspan=3, columnspan=2, padx=5, pady=(5,5), sticky='nsew')
    search_frame.grid_propagate(False)

    subframe1 = Frame(search_frame)
    subframe1.grid(row=0, column=0, pady=(20, 10), sticky='nsew')
    subframe1.grid_propagate(False)
    mass_label = Label(subframe1, text="<= peptide mass <=")
    mass_label.grid(row=0, column=2, padx=10, sticky='new')
    lower_en = Entry(subframe1, width=5)
    lower_en.insert(0, 600)
    lower_en.bind("<FocusIn>", lambda event: handle_focus_in(key=lower_en))
    lower_en.grid(row=0, column=0, padx=(10, 0), sticky='new')

    higher_en = Entry(subframe1, width=5)
    higher_en.insert(0, 6000)
    higher_en.bind("<FocusIn>", lambda event: handle_focus_in(key=higher_en))
    higher_en.grid(row=0, column=3, sticky='new')

    leftda_label = Label(subframe1, text="Da")
    leftda_label.grid(row=0, column=1, sticky='new')
    rightda_label = Label(subframe1, text="Da")
    rightda_label.grid(row=0, column=4, sticky='new')

    parse_label = Label(subframe1, text="Digestion rule :")
    parse_en = ttk.Combobox(subframe1, width=8, state="readonly")
    parse_en['values'] = ("trypsin")
    parse_en.current(0)
    parse_label.grid(row=0, column=5, padx=(50, 5), sticky='new')
    parse_en.grid(row=0, column=6, padx=(0,10), sticky='new')
    subframe1.grid_columnconfigure(0, weight=1)
    subframe1.grid_columnconfigure(1, weight=1)
    subframe1.grid_columnconfigure(2, weight=1)
    subframe1.grid_columnconfigure(3, weight=1)
    subframe1.grid_columnconfigure(4, weight=1)
    subframe1.grid_columnconfigure(5, weight=1)
    subframe1.grid_columnconfigure(6, weight=1)

    subframe2 = Frame(search_frame)
    subframe2.grid(row=1, column=0, pady=10, sticky='nsew')
    subframe2.grid_propagate(False)
    miss_label = Label(subframe2, text="Miss cleavage :")
    miss_en = ttk.Combobox(subframe2, width=2, state="readonly")
    miss_en['values'] = (0, 1, 2, 3)
    miss_en.current(2)
    miss_label.grid(row=0, column=0, padx=(5, 5), sticky='new')
    miss_en.grid(row=0,column=1, sticky='new')

    ms1tol_label = Label(subframe2, text="MS1 tolerance :")
    ms1tol_en = Entry(subframe2, width=4)
    ms1tol_en.insert(0, 20)
    ms1unit_label = Label(subframe2, text="ppm")
    ms1tol_label.grid(row=0, column=2, padx=(10, 5), sticky='new')
    ms1tol_en.grid(row=0,column=3, sticky='new')
    ms1tol_en.bind("<FocusIn>", lambda event: handle_focus_in(key=ms1tol_en))
    ms1unit_label.grid(row=0, column=4, sticky='new')

    ms2tol_label = Label(subframe2, text="MS2 tolerance :")
    ms2tol_en = Entry(subframe2, width=4)
    ms2tol_en.insert(0, 0.02)
    ms2unit_label = Label(subframe2, text="Da")
    ms2tol_label.grid(row=0, column=5, padx=(10, 5), sticky='new')
    ms2tol_en.grid(row=0,column=6, sticky='new')
    ms2tol_en.bind("<FocusIn>", lambda event: handle_focus_in(key=ms2tol_en))
    ms2unit_label.grid(row=0, column=7, padx=(0,10), sticky='new')

    subframe2.grid_columnconfigure(0, weight=1)
    subframe2.grid_columnconfigure(1, weight=1)
    subframe2.grid_columnconfigure(2, weight=1)
    subframe2.grid_columnconfigure(3, weight=1)
    subframe2.grid_columnconfigure(4, weight=1)
    subframe2.grid_columnconfigure(5, weight=1)
    subframe2.grid_columnconfigure(6, weight=1)
    subframe2.grid_columnconfigure(7, weight=1)

    subframe3 = Frame(search_frame)
    subframe3.grid(row=2, column=0, pady=10, sticky='nsew')
    fdr_label = Label(subframe3, text="FDR setting :")
    fdr_en = Entry(subframe3, width=2)
    fdr_en.insert(0, 1)
    fdr_en.bind("<FocusIn>", lambda event: handle_focus_in(key=fdr_en))

    fdrunit_label = Label(subframe3, text="%")
    fdr_label.grid(row=0, column=0, padx=5)
    fdr_en.grid(row=0, column=1)
    fdrunit_label.grid(row=0,column=2)

    search_frame.grid_columnconfigure(0, weight=1)
    search_frame.grid_rowconfigure(0, weight=1)
    search_frame.grid_rowconfigure(1, weight=1)
    search_frame.grid_rowconfigure(2, weight=1)

    '''linker frame'''
    linkerDict = {}
    linker_frame = LabelFrame(non_win, text='Cross-linker')
    linker_frame.grid(row=0, column=2, rowspan=2, columnspan=2, padx=5, pady=(10,5), sticky='nsew')
    linker_frame.grid_propagate(False)
    nonlinker_label = Label(linker_frame, text='Select non-cleavable cross linker :')
    nonlinker_en = ttk.Combobox(linker_frame, width=4, state="readonly")
    refresh_nonlinker()
    nonlinker_label.grid(row=0, column=0, columnspan=3, padx=5, pady=5, sticky='nsew')
    nonlinker_en.grid(row=0, column=3, padx=20, pady=10, sticky='nsew')
    db_frame = LabelFrame(linker_frame, text='Non-cleavable cross linker database')
    db_frame.grid(row=1, column=0, columnspan=3, padx=5, pady=(10, 5), sticky='nsew')
    db_frame.grid_propagate(False)
    add_btn = Button(db_frame, text='Add', width=8, command=add_nonlinker)
    edit_btn = Button(db_frame, text='Edit', width=8, command=edit_nonlinker)
    del_btn = Button(db_frame, text='Delete', width=8, command=delete_nonlinker)
    ref_btn = Button(linker_frame, text='Update database', width=10, command=refresh_nonlinker)
    add_btn.grid(row=0, column=0, padx=10, pady=5, sticky='nsew')
    add_btn.bind("<Enter>", on_enter_add_btn)
    add_btn.bind("<Leave>", on_non_leave)
    edit_btn.grid(row=0, column=1, padx=10, pady=5, sticky='nsew')
    edit_btn.bind("<Enter>", on_enter_edit_btn)
    edit_btn.bind("<Leave>", on_non_leave)
    del_btn.grid(row=0, column=2, padx=10, pady=5, sticky='nsew')
    del_btn.bind("<Enter>", on_enter_del_btn)
    del_btn.bind("<Leave>", on_non_leave)
    ref_btn.grid(row=1, column=3, padx=20, pady=(20, 10), sticky='nsew')
    ref_btn.bind("<Enter>", on_enter_ref_btn)
    ref_btn.bind("<Leave>", on_non_leave)

    db_frame.grid_columnconfigure(0, weight=1)
    db_frame.grid_columnconfigure(1, weight=1)
    db_frame.grid_columnconfigure(2, weight=1)
    db_frame.grid_rowconfigure(0, weight=1)

    linker_frame.grid_rowconfigure(0, weight=1)
    linker_frame.grid_rowconfigure(1, weight=1)
    linker_frame.grid_columnconfigure(0, weight=1)
    linker_frame.grid_columnconfigure(1, weight=1)
    linker_frame.grid_columnconfigure(2, weight=1)
    linker_frame.grid_columnconfigure(3, weight=1)

    '''mod frame'''
    mod_frame = LabelFrame(non_win, text='Modifications')
    mod_frame.grid(row=2, column=2, rowspan=4, columnspan=2, padx=5, pady=(5,5), sticky='nsew')
    mod_frame.grid_propagate(False)
    mod_db_frame = LabelFrame(mod_frame, text='Mods. database')
    mod_db_frame.grid(row=0, column=3, rowspan=3, padx=5, pady=(10, 5), sticky='nsew')
    mod_db_frame.grid_propagate(False)

    fixed_label = Label(mod_frame, text="fixed mods.")
    fixed_en = Listbox(mod_frame, height=7, width=15)
    var_label = Label(mod_frame, text="variable mods.")
    var_en = Listbox(mod_frame, height=7, width=15)
    fixed_label.grid(row=0, column=0, pady=5, sticky='nsew')
    fixed_en.grid(row=1, column=0, padx=5, pady=5, sticky='nsew')
    var_label.grid(row=2, column=0, pady=5, sticky='nsew')
    var_en.grid(row=3, column=0, padx=5, pady=5, sticky='nsew')

    add_fixed_btn = Button(mod_frame, text="<", command=addfixedmod, anchor='n')
    rem_fixed_btn = Button(mod_frame, text=">", command=remfixedmod, anchor='s')
    add_var_btn = Button(mod_frame, text="<", command=addvarmod, anchor='n')
    rem_var_btn = Button(mod_frame, text=">", command=remvarmod, anchor='s')
    add_fixed_btn.grid(row=1, column=1, pady=15, sticky='n')
    add_fixed_btn.bind("<Enter>", on_enter_add_fixed_btn)
    add_fixed_btn.bind("<Leave>", on_non_leave)
    rem_fixed_btn.grid(row=1, column=1, pady=15, sticky='s')
    rem_fixed_btn.bind("<Enter>", on_enter_rem_fixed_btn)
    rem_fixed_btn.bind("<Leave>", on_non_leave)
    add_var_btn.grid(row=3, column=1, pady=15, sticky='n')
    add_var_btn.bind("<Enter>", on_enter_add_var_btn)
    add_var_btn.bind("<Leave>", on_non_leave)
    rem_var_btn.grid(row=3, column=1, pady=15, sticky='s')
    rem_var_btn.bind("<Enter>", on_enter_rem_var_btn)
    rem_var_btn.bind("<Leave>", on_non_leave)

    modDict = {}

    mod_en = Listbox(mod_frame, height=21, width=15)
    refresh_mod_nonlinker()
    mod_en.grid(row=0, column=2, rowspan=4, padx=5, pady=5, sticky='nsew')

    mod_add_btn = Button(mod_db_frame, text='Add', width=8, command=add_mod_nonlinker)
    mod_edit_btn = Button(mod_db_frame, text='Edit', width=8, command=edit_mod_nonlinker)
    mod_del_btn = Button(mod_db_frame, text='Delete', width=8, command=delete_mod_nonlinker)
    mod_ref_btn = Button(mod_frame, text='Update mods.', width=10, command=refresh_mod_nonlinker)
    mod_add_btn.grid(row=0, column=0, padx=5, pady=5, sticky='nsew')
    mod_add_btn.bind("<Enter>", on_enter_mod_add_btn)
    mod_add_btn.bind("<Leave>", on_non_leave)
    mod_edit_btn.grid(row=1, column=0, padx=5, pady=5, sticky='nsew')
    mod_edit_btn.bind("<Enter>", on_enter_mod_edit_btn)
    mod_edit_btn.bind("<Leave>", on_non_leave)
    mod_del_btn.grid(row=2, column=0, padx=5, pady=5, sticky='nsew')
    mod_del_btn.bind("<Enter>", on_enter_mod_del_btn)
    mod_del_btn.bind("<Leave>", on_non_leave)
    mod_ref_btn.grid(row=3, column=3, padx=5, pady=(25, 25), sticky='nsew')
    mod_ref_btn.bind("<Enter>", on_enter_mod_ref_btn)
    mod_ref_btn.bind("<Leave>", on_non_leave)

    mod_db_frame.grid_columnconfigure(0, weight=1)
    mod_db_frame.grid_rowconfigure(0, weight=1)
    mod_db_frame.grid_rowconfigure(1, weight=1)
    mod_db_frame.grid_rowconfigure(2, weight=1)

    mod_frame.grid_rowconfigure(0, weight=1)
    mod_frame.grid_rowconfigure(1, weight=3)
    mod_frame.grid_rowconfigure(2, weight=1)
    mod_frame.grid_rowconfigure(3, weight=3)
    mod_frame.grid_columnconfigure(0, weight=2)
    mod_frame.grid_columnconfigure(1, weight=1)
    mod_frame.grid_columnconfigure(2, weight=3)
    mod_frame.grid_columnconfigure(3, weight=2)

    '''back, stop, start'''

    img_back = ImageTk.PhotoImage(Image.open(r"arrow.png").resize((25, 20)))
    back_btn = Button(non_win, text="Back", image=img_back, borderwidth=2, compound=LEFT, command=back2normal1)
    back_btn.grid(row=6, column=0, padx=5, ipadx=10, ipady=5, sticky='w')
    back_btn.bind("<Enter>", on_enter_back_btn)
    back_btn.bind("<Leave>", on_non_leave)
    img_stop = ImageTk.PhotoImage(Image.open(r"stop.png").resize((20, 20)))
    stop_btn = Button(non_win, text='Stop', image=img_stop, borderwidth=2, compound=LEFT, state=DISABLED, command=stop)
    stop_btn.grid(row=6, column=3, padx=5, ipadx=10, ipady=5, sticky='w')
    stop_btn.bind("<Enter>", on_enter_stop_btn)
    stop_btn.bind("<Leave>", on_non_leave)
    img_start = ImageTk.PhotoImage(Image.open(r"start.png").resize((25, 20)))
    start_btn = Button(non_win, text='Start', image=img_start, borderwidth=2, compound=LEFT, command=run_xolik)
    start_btn.grid(row=6, column=3, padx=5, ipadx=10, ipady=5, sticky='e')
    start_btn.bind("<Enter>", on_enter_start_btn)
    start_btn.bind("<Leave>", on_non_leave)

    '''print frame'''
    print_frame = Frame(non_win)
    print_frame.grid(row=7, column=0, rowspan=3, columnspan=4, padx=5, pady=(5,5), sticky='nsew')
    print_frame.grid_propagate(False)
    print_en = Text(print_frame, bg="white", relief="sunken", wrap=NONE, state=DISABLED)
    print_en.grid(row=0, column=0, sticky='nsew')

    print_frame.grid_columnconfigure(0, weight=1)
    print_frame.grid_rowconfigure(0, weight=1)
    '''status'''
    non_status = Label(non_win, text='Non-cleavable cross-linker searching module', bd=1, relief=SUNKEN, anchor=E)
    non_status.grid(row=10, column=0, columnspan=4, ipady=4, sticky='sew')

    non_win.grid_columnconfigure(0, weight=1)
    non_win.grid_columnconfigure(1, weight=1)
    non_win.grid_columnconfigure(2, weight=1)
    non_win.grid_columnconfigure(3, weight=1)

    non_win.grid_rowconfigure(0, weight=1)
    non_win.grid_rowconfigure(1, weight=1)
    non_win.grid_rowconfigure(2, weight=1)
    non_win.grid_rowconfigure(3, weight=1)
    non_win.grid_rowconfigure(4, weight=1)
    non_win.grid_rowconfigure(5, weight=1)
    non_win.grid_rowconfigure(7, weight=1)
    non_win.grid_rowconfigure(8, weight=1)
    non_win.grid_rowconfigure(9, weight=1)


def on_non_leave(e):
    """Noncleaevable module: when leave button"""
    e.widget['background'] = 'SystemButtonFace'
    non_status.config(text='Non-cleavable cross-linker searching module')


def on_non_add_linker_leave(e):
    """Noncleavable: add linker frame"""
    add_status.config(text='Add non-cleavable cross linker')


def on_non_edit_linker_leave(e):
    """Noncleavable: edit linker frame"""
    add_status.config(text='Edit non-cleavable cross linker')


def on_non_add_mod_leave(e):
    """Noncleavable module: enter add modification button"""
    add_status.config(text='Add modification')


def on_non_edit_mod_leave(e):
    """Noncleavable module: leave modification frame"""
    add_status.config(text='Edit modification')

'''##################################################### Root #######################################################'''


root = Tk()
root.title("ECL")
root.iconbitmap(r"ECL.ico")
root_w = 900
root_h = 600
screen_width = root.winfo_screenwidth()
screen_height = root.winfo_screenheight()
root_x = int(screen_width/2 - root_w/2)
root_y = int(screen_height/2 - root_h/2)
root.geometry('%dx%d+%d+%d' % (root_w, root_h, root_x, root_y))
# root.columnconfigure(1, minsize=100, weight=0)
fasta = ""
cleav_in = 0
file_in = 1
'''create database'''
# conn = sqlite3.connect('paramsDB.db')
# c = conn.cursor()
# c.execute("""CREATE TABLE nonlinker (
#         linker text,
#         site text,
#         mass real
#         )""")
#
# c.execute("""CREATE TABLE nonlinkermod (
#         name text,
#         site text,
#         mass real
#         )""")
#
# c.execute("""CREATE TABLE cleavelinker (
#         linker text,
#         site text,
#         mass real,
#         m_short real,
#         m_long real
#         )""")
#
# c.execute("""CREATE TABLE cleavelinkermod (
#         name text,
#         site text,
#         mass real
#         )""")
#
# conn.commit()
# conn.close()
'''root'''
title = Label(root, text='Welcome to ECL 3.0!', font=('Helvetica', 40))
title.grid(row=0, column=0, columnspan=4)

my_img1 = ImageTk.PhotoImage(Image.open(r"ECL.png").resize((425, 200)))

my_icon = Label(image=my_img1)
my_icon.grid(row=1, column=0, rowspan=5, columnspan=2)
intro = Label(root, text='Introduction:', font=('Helvetica', 20, 'bold'))
intro.grid(row=1, column=2, columnspan=2, sticky='w')

intro_text = 'ECL is a mass spectrometry search engine developed from Yu\'s group at HKUST. ' \
             'It aims to provide the solution to the cross-linking mass spectrometry (XL-MS) data analysis. ' \
             'It exhaustively searches all the possible cross-links in a linear time complexity and ' \
             'integrates the protein feedback mechanism to improve the software sensitivity.'

intro_content = Label(root, text=intro_text,  wraplengt=420, justify=LEFT)
intro_content.grid(row=2, column=2, columnspan=2)

cite = Label(root, text='Cite us:', font=('Helvetica', 20, 'bold'))
cite.grid(row=3, column=2, columnspan=2, sticky='w')

cite_text1 = 'Chen Zhou, et al. "Exhaustive Cross-Linking Search with Protein Feedback", Journal of Proteome Research, 2022.'
cite_content1 = Label(root, text=cite_text1, wraplengt=400, justify=LEFT, cursor='hand2')
cite_content1.grid(row=4, column=2, columnspan=2, pady=5)
cite_content1.bind("<Button-1>", lambda e: callurl("https://doi.org/10.1021/acs.jproteome.2c00500"))

cite_text2 = 'Jiaan Dai, et al. "Xolik: Finding Cross-Linked Peptides with Maximum Paired Scores in Linear Time", Bioinformatics, 2019.'
cite_content2 = Label(root, text=cite_text2, wraplengt=400, justify=LEFT, cursor='hand2')
cite_content2.grid(row=5, column=2, columnspan=2, pady=5)
cite_content2.bind("<Button-1>", lambda e: callurl("https://doi.org/10.1093/bioinformatics/bty526"))

cite_text3 = 'Fengchao Yu, et al. "Exhaustively Identifying Cross-Linked Peptides with a Linear Computational Complexity", Journal of Proteome Research, 2017.'
cite_content3 = Label(root, text=cite_text3, wraplengt=400, justify=LEFT, cursor='hand2')
cite_content3.grid(row=6, column=2, columnspan=2, pady=5)
cite_content3.bind("<Button-1>", lambda e: callurl("https://doi.org/10.1021/acs.jproteome.7b00338"))

cite_text4 = 'Fengchao Yu, et al. "ECL: An Exhaustive Search Tool for The Identification of Cross-Linked Peptides Using Whole Database", BMC Bioinformatics, 2016.'
cite_content4 = Label(root, text=cite_text4, wraplengt=400, justify=LEFT, cursor='hand2')
cite_content4.grid(row=7, column=2, columnspan=2, pady=5)
cite_content4.bind("<Button-1>", lambda e: callurl("https://doi.org/10.1186/s12859-016-1073-y"))

'''contact info'''
contact_info = Label(root, text='Contact Info: czhouau@connect.ust.hk', font=('Arial', 9))
contact_info.grid(row=8, column=2, columnspan=2, sticky='w')

'''cleavable button in root'''
cleavable = Button(root, text = "Cleavable\ndata analysis", width=10, font=('Helvetica', 15), borderwidth=5, relief=RAISED, command=cleav)
cleavable.grid(row=6, column=1, padx=20, pady=30, rowspan=3, sticky="nsew")
cleavable.bind("<Enter>", on_enter_cleavable)
cleavable.bind("<Leave>", on_leave)

'''non cleavable button in root'''
noncleavable = Button(root, text = "Non-cleavable\ndata analysis", width=10, font=('Helvetica', 15), borderwidth=5, relief=RAISED, command=noncleav)
noncleavable.grid(row=6, column=0, padx=20, pady=30, rowspan=3, sticky="nsew")
noncleavable.bind("<Enter>", on_enter_noncleavable)
noncleavable.bind("<Leave>", on_leave)


'''scalable button'''
root.grid_columnconfigure(0, weight=1)
root.grid_columnconfigure(1, weight=1)
root.grid_columnconfigure(2, weight=1)
root.grid_columnconfigure(3, weight=1)

root.grid_rowconfigure(0, weight=2)
root.grid_rowconfigure(1, weight=1)
root.grid_rowconfigure(2, weight=2)
root.grid_rowconfigure(3, weight=1)
root.grid_rowconfigure(8, weight=1)
root.grid_rowconfigure(9, weight=1)


status_label = Label(root, text='Welcome to ECL 3.0!', bd=1, relief=SUNKEN, anchor=E)
status_label.grid(row=9, column=0, columnspan=4, ipady=4, sticky='sew')

root.mainloop()
