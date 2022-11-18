#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
"""
mzcalc
20 october 2022
"""
#
import wx
import os
from main_aquasearch import main_aquasearch
import score_matchms as sm
import PCA_analysis 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
import protein_signals as ps
import pd_maldi_match as pdmm
import standard_incorporation as si
import sqlite_requests as sr
import pandas as pd 
import sqlite3
import db_elimination as dbe
import help_strings as hs



#
class AquasearchFrame(main_aquasearch):
    
    def __init__(self, *args, **kwds):

        main_aquasearch.__init__(self, *args, **kwds)
        
        self.dic_s = {}
        self.Bind(wx.EVT_BUTTON, self.on_bt_query, self.button_run_query)
        self.Bind(wx.EVT_BUTTON, self.on_bt_plot, self.button_pca)
        self.Bind(wx.EVT_BUTTON, self.on_bt_1pt_mix, self.btn_1lbl_mix_run)
        self.Bind(wx.EVT_BUTTON, self.on_bt_1pt_std, self.btn_1lbl_std_run)
        self.Bind(wx.EVT_BUTTON, self.on_bt_2pt_mix, self.btn_2lbl_mix_run)
        self.Bind(wx.EVT_BUTTON, self.on_bt_2pt_std, self.btn_2lbl_std_run) 
        self.Bind(wx.EVT_BUTTON, self.DownPlot, self.button_download)
        self.Bind(wx.EVT_BUTTON, self.on_bt_del_prot, self.btn_del_prot)
        self.Bind(wx.EVT_BUTTON, self.on_bt_del_sam, self.btn_del_sam)
        
        self.Bind(wx.EVT_COMBOBOX, self.OnSelectionProtein, self.cb_database_delete_prot)
        self.Bind(wx.EVT_COMBOBOX, self.OnSelectionSample, self.cb_database_delete_sam)
        self.Bind(wx.EVT_COMBOBOX, self.NewDb, self.cb_database_ne)
        
        self.Bind(wx.EVT_SPINCTRLDOUBLE, self.ComponentsSel, self.sc_number_PCs)
        
        self.Bind(wx.EVT_SIZE,self.OnSize)
        self.Bind(wx.EVT_IDLE,self.OnIdle)
        
        self.tolerance_text.SetToolTip(hs.tolerance)
        self.database_text.SetToolTip(hs.database)
        self.advanced_text.SetToolTip(hs.advanced_details)
        self.database_txt_ne.SetToolTip(hs.database_ne)
        
    def NewDb(self, evt):
        try:
            db = self.cb_database_ne.GetSelection()
            db = self.dbs_[db]
            if db == "New database":
                new_db = ''
                while not new_db:
                    dlg = wx.TextEntryDialog(self,'New database', 'Database name', "", wx.OK|wx.CANCEL)
                    answ = dlg.ShowModal()
                    if answ == wx.ID_OK:
                        new_db = dlg.GetValue()
                        if not new_db:
                            wx.MessageBox('Error: Introduce the name of the new database', 'Info', wx.OK | wx.ICON_INFORMATION)
                    
                    elif answ == wx.ID_CANCEL:  
                        break
                        
                if not new_db:
                    dlg.Close()
                else:
                    sr.create_db_by_user(new_db)
                    new_list = [new_db] + self.dbs_
                
                    self.cb_database_ne.Clear()
                    for l in new_list:
                        self.cb_database_ne.Append(l)
                        
                    new_list.remove("New database")
                    self.cb_database_delete_prot.Clear()
                    for l in new_list:
                        self.cb_database_delete_prot.Append(l)
                        
                    self.cb_database_delete_sam.Clear()
                    for l in new_list:
                        self.cb_database_delete_sam.Append(l)
                
                    self.cb_database.Clear()
                    for l in new_list:
                        self.cb_database.Append(l)
                
        except sqlite3.OperationalError: 
            wx.MessageBox('Error: No database selected', 'Info', wx.OK | wx.ICON_INFORMATION)
                
    
    def on_bt_del_prot(self, evt):
        try:
            db_ = self.cb_database_delete_prot.GetValue()
            protein = self.cb_delete_prot.GetValue()
            if db_ and protein:
            
                dlg = wx.MessageDialog(None, 'Are you sure you want to delete the protein "' + protein + '"?', 'Confirm Delete',
                                       wx.YES_NO | wx.ICON_QUESTION)
                dlg.SetYesNoLabels('Yes', 'No')
                answ = dlg.ShowModal()
                if answ == wx.ID_YES:
                    wx.BeginBusyCursor(wx.Cursor(wx.CURSOR_ARROWWAIT))
                    dbe.protein_elimination(db_, protein)
                    self.OnSelectionProtein(self.cb_database_delete_prot)
                    wx.EndBusyCursor()
                    wx.MessageBox('Process has been successfully completed !!', 'Info', wx.OK | wx.ICON_INFORMATION)
                elif answ == wx.ID_NO:
                    dlg.Close()
            else:
                dbe.protein_elimination(db_, protein)
                
        except sqlite3.OperationalError: 
            wx.MessageBox('Error: No database selected', 'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
        except IndexError:
            wx.MessageBox('Error: No protein selected', 'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
        
    def on_bt_del_sam(self, evt):
        try:
            db_ = self.cb_database_delete_sam.GetValue()
            sample = self.cb_delete_sam.GetValue()
            if db_ and sample:             
                dlg = wx.MessageDialog(None, 'Are you sure you want to delete the sample "'+ sample + '"?',
                                       'Confirm Delete', wx.YES_NO | wx.ICON_QUESTION)
                dlg.SetYesNoLabels('Yes', 'No')
                answ = dlg.ShowModal()
                if answ == wx.ID_YES:
                    wx.BeginBusyCursor(wx.Cursor(wx.CURSOR_ARROWWAIT))
                    dbe.sample_elimination(db_, sample)
                    self.OnSelectionSample(self.cb_database_delete_sam)
                    wx.EndBusyCursor()
                    wx.MessageBox('Process has been successfully completed !!', 'Info', wx.OK | wx.ICON_INFORMATION)
                elif answ == wx.ID_NO:
                    dlg.Close()
            else:
                if not db_:
                    dbe.sample_elimination(db_, sample)
                else:
                    raise IndexError
                    
        
        except sqlite3.OperationalError: 
            wx.MessageBox('Error: No database selected', 'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
        except IndexError:
            wx.MessageBox('Error: No sample selected', 'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
    
    def OnSelectionProtein(self, evt):
        db_chosen = self.cb_database_delete_prot.GetSelection()
        db_chosen = self.dbs_option[db_chosen]
        if 'Aquasearch_study' not in db_chosen:
            db_chosen = 'Aquasearch_study - ' + db_chosen
            
        table_ = pd.DataFrame(sr.table_download(db_chosen, 'Protein_codes'))
        
        self.cb_delete_prot.Clear()
        for cod in table_.iloc[:,1]:
            self.cb_delete_prot.Append(cod)
    
    def OnSelectionSample(self, evt):
        db_chosen = self.cb_database_delete_sam.GetSelection()
        db_chosen = self.dbs_option[db_chosen]
        if 'Aquasearch_study' not in db_chosen:
            db_chosen = 'Aquasearch_study - ' + db_chosen
            
        table_ = pd.DataFrame(sr.table_download(db_chosen, 'Spectrums_table'))
        table_ = table_.iloc[:,4].unique()
        table_.sort()
        
        self.cb_delete_sam.Clear()
        for cod in table_:
            self.cb_delete_sam.Append(cod)
    
    def ComponentsSel(self, evt):
        n_pc = int(self.sc_number_PCs.GetValue())
        
        self.PC_x.Clear()
        self.PC_y.Clear()
        self.PC_z.Clear()
        
        self.PC_z.Append('None')
        for n in range(n_pc):
            self.PC_x.Append(str(n+1))
            self.PC_y.Append(str(n+1))
            self.PC_z.Append(str(n+1))
    
    def OnSize(self, evt):
        self.resized = True
        evt.Skip()
        
    def OnIdle(self, evt):
        if self.resized: 
            self.CleanPlot()
            self.figure.set_size_inches(self.panel_canvas.Size[0]/105, self.panel_canvas.Size[1]/105)
            self.canvas = FigureCanvas(self.panel_canvas, -1,  self.figure)
            
            self.resized = False
    
    def on_bt_query(self, evt):
        try:
            wx.BeginBusyCursor(wx.Cursor(wx.CURSOR_ARROWWAIT))
            self.CleanResults()
        
            self.dic_s = self.Calculate_dict()
        
            for numero in self.dic_s.keys():
                self.lc_samples.InsertItem(1,numero)
            
            self.Bind(wx.EVT_LIST_ITEM_ACTIVATED, lambda evt, temp=self.dic_s: self.results_display(evt, temp))
            wx.EndBusyCursor()
            
        except FileNotFoundError:
            wx.MessageBox('Error: No file pathway selected', 'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
        except sqlite3.OperationalError: 
            wx.MessageBox('Error: No database selected', 'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
        
    def on_bt_plot(self, evt):
        try: 
            self.CleanPlot()
        
            if not self.dic_s:
                self.dic_s = self.Calculate_dict()
        
            n_pcs = int(self.sc_number_PCs.GetValue())
            pc_x = int(self.PC_x.GetValue())
            pc_y = int(self.PC_y.GetValue())
            pc_z = self.PC_z.GetValue()
            # title = str(self.canvas.GetSize())
            title = self.fig_title.GetValue()
            point_size = self.point_size.GetValue()
        
            if pc_z == 'None':
                data_matrix, vari, samp_nam = PCA_analysis.PCA_2d(self.dic_s, n_pcs)
                self.draw_2d(self, data_matrix, vari, pc_x, pc_y, title, point_size, samp_nam) 
            
            
            else:
                data_matrix, vari, samp_nam = PCA_analysis.PCA_3d(self.dic_s, n_pcs)
                self.draw_3d(self, data_matrix, vari, pc_x, pc_y, int(pc_z), title, point_size, samp_nam)
                
        except sqlite3.OperationalError: 
            wx.MessageBox('Error: No database selected', 'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
        except FileNotFoundError:
            wx.MessageBox('Error: No file pathway selected', 'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
        except ValueError:
            wx.MessageBox('Error: Select all PC correctly', 'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
            
    def on_bt_1pt_mix(self, evt):
        try:
            wx.BeginBusyCursor(wx.Cursor(wx.CURSOR_ARROWWAIT))
            db = self.cb_database_ne.GetValue()
            if db:
                if not 'Aquasearch_study' in db:
                    db = 'Aquasearch_study - ' + db
        
            maldi = self.btn_1lbl_mix_maldi.GetValue()
            peptide = self.btn_1lbl_mix_pep.GetValue()
            protein = self.btn_1lbl_mix_prot.GetValue()
        
            ppm = self.sc_1lbl_mix_ppm.GetValue()
        
            uni_pep = self.cb_1prot_uniq.GetValue()   
            if uni_pep == "Selection ":
                uni_pep = 0
            elif uni_pep == "All posibilities":
                uni_pep = 1
            else:
                raise ValueError
            
            code = self.txt_1prot_mix_code.GetValue()
            rank = self.sc_1prot_rank.GetValue()
            sample_name = self.txt_1prot_mix_name.GetValue()
            if not code or not sample_name:
                raise ValueError
        
            test_pdmm = pdmm.txt_complete(maldi, peptide, protein, rank, ppm, uni_pep)
        
            ps.fill_table(code, test_pdmm, sample_name, db)
            
            wx.EndBusyCursor()
            wx.MessageBox('Process has been successfully completed !!', 'Info', wx.OK | wx.ICON_INFORMATION)
        except sqlite3.OperationalError: 
            wx.MessageBox('Error: No database selected', 'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
        except FileNotFoundError:
            wx.MessageBox('Error: No file pathway selected', 'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
        except ValueError:
            wx.MessageBox('Error: Select all the parameters', 'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
        
    def on_bt_1pt_std(self, evt):
        try:
            wx.BeginBusyCursor(wx.Cursor(wx.CURSOR_ARROWWAIT))
            db = self.cb_database_ne.GetValue()
            if db:
                if not 'Aquasearch_study' in db:
                    db = 'Aquasearch_study - ' + db
        
            maldi = self.btn_1lbl_std_maldi.GetValue()
            peptide = self.btn_1lbl_std_uni.GetValue()
        
            ppm = self.sc_1lbl_std_ppm.GetValue()
            code = self.txt_1prot_std_code.GetValue()
            sample_name = self.txt_1prot_std_name.GetValue()
            if not code or not sample_name:
                raise ValueError
        
            si.standard_complete(maldi, peptide, code, sample_name, ppm, db)
            
            wx.EndBusyCursor()
            wx.MessageBox('Process has been successfully completed !!', 'Info', wx.OK | wx.ICON_INFORMATION)
        except sqlite3.OperationalError: 
            wx.MessageBox('Error: No database selected', 'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
        except FileNotFoundError:
            wx.MessageBox('Error: No file pathway selected', 'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
        except ValueError:
            wx.MessageBox('Error: Select all the parameters', 'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
        
    def on_bt_2pt_mix(self, evt):
        try:
            wx.BeginBusyCursor(wx.Cursor(wx.CURSOR_ARROWWAIT))
            db = self.cb_database_ne.GetValue()
            if db:
                if not 'Aquasearch_study' in db:
                    db = 'Aquasearch_study - ' + db
        
            maldi = self.btn_2lbl_mix_maldi.GetValue()
            peptide = self.btn_2lbl_mix_pep.GetValue()
            protein = self.btn_2lbl_mix_prot.GetValue()
        
            ppm = self.sc_2lbl_mix_ppm.GetValue()
            uni_pep = self.cb_2prot_uniq.GetValue()
            if uni_pep == "Selection ":
                uni_pep = 0
            elif uni_pep == "All posibilities":
                uni_pep = 1
            else:
                raise ValueError
            
            code1 = self.txt_2prot_mix_code1.GetValue()
            code2 = self.txt_2prot_mix_code2.GetValue()
            new_entry = self.txt_2prot_mix_new.GetValue()
            rank = self.sc_2prot_rank.GetValue()
            sample_name = self.txt_2prot_mix_name.GetValue()
            if not code1 or not code2 or not new_entry or not sample_name:
                raise ValueError
        
            test_pdmm = pdmm.txt_complete(maldi, peptide, protein, rank, ppm, uni_pep)
            ps.table_union(new_entry, code1,code2, test_pdmm, sample_name, db)
            
            wx.EndBusyCursor()
            wx.MessageBox('Process has been successfully completed !!', 'Info', wx.OK | wx.ICON_INFORMATION)
        except sqlite3.OperationalError: 
            wx.MessageBox('Error: No database selected', 'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
        except FileNotFoundError:
            wx.MessageBox('Error: No file pathway selected', 'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
        except ValueError:
            wx.MessageBox('Error: Select all the parameters', 'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
    
    def on_bt_2pt_std(self, evt):
        try:
            wx.BeginBusyCursor(wx.Cursor(wx.CURSOR_ARROWWAIT))
            db = self.cb_database_ne.GetValue()
            if db:
                if not 'Aquasearch_study' in db:
                    db = 'Aquasearch_study - ' + db
        
            maldi = self.btn_2lbl_std_maldi.GetValue()
            peptide = self.btn_2lbl_std_uni.GetValue()
        
            ppm = self.sc_2lbl_std_ppm.GetValue()
            code = self.txt_2prot_std_code.GetValue()
            new_entry = self.txt_2prot_std_new.GetValue()
            sample_name = self.txt_2prot_std_name.GetValue()
            if not code or not new_entry or not sample_name:
                raise ValueError
        
            si.standard_complete_union(maldi, peptide, code, new_entry, sample_name, ppm, db)
            
            wx.EndBusyCursor()
            wx.MessageBox('Process has been successfully completed !!', 'Info', wx.OK | wx.ICON_INFORMATION)
            
        except sqlite3.OperationalError: 
            wx.MessageBox('Error: No database selected', 'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
        except FileNotFoundError:
            wx.MessageBox('Error: No file pathway selected', 'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
        except ValueError:
            wx.MessageBox('Error: Select all the parameters', 'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
            
    def draw_3d(self, evt, data_matrix, var, x, y, z, title, size_point, sample_names):
    
        self.figure = plt.figure(figsize=(self.panel_canvas.Size[0]/105, self.panel_canvas.Size[1]/105))
        # self.figure = plt.figure()
        ax = self.figure.add_subplot(projection='3d')
        ax.scatter(data_matrix[:,x-1], data_matrix[:,y-1], data_matrix[:,z-1], s = size_point)
        
        for i, txt in enumerate(sample_names): 
            ax.text3D(data_matrix[i,x-1]+0.015, data_matrix[i,y-1]+0.015, data_matrix[i,z-1]+0.015,txt)
            
        ax.set_xlabel('PC' + str(x) + ' (' + str(round(var[x-1]*100, 2)) + '%)')
        ax.set_ylabel('PC' + str(y) + ' (' + str(round(var[y-1]*100, 2)) + '%)')
        ax.set_zlabel('PC' + str(z) + ' (' + str(round(var[z-1]*100, 2)) + '%)')
        plt.title(title)
        self.canvas = FigureCanvas(self.panel_canvas, -1, self.figure)
    
    def draw_2d(self, evt, data_matrix, var, x, y, title, size_point, sample_names):
        self.figure = plt.figure(figsize=(self.panel_canvas.Size[0]/105, self.panel_canvas.Size[1]/105))
        # self.figure = plt.figure()
        ax = self.figure.add_subplot()
        ax.scatter(data_matrix[:,x-1], data_matrix[:,y-1], s = size_point)
        
        for i, txt in enumerate(sample_names):
            ax.annotate(txt, (data_matrix[i,x-1]+0.015, data_matrix[i,y-1]+0.015))
            
        ax.set_xlabel('PC' + str(x) + ' (' + str(round(var[x-1]*100, 2)) + '%)')
        ax.set_ylabel('PC' + str(y) + ' (' + str(round(var[y-1]*100, 2)) + '%)')
        plt.title(title)
        self.canvas = FigureCanvas(self.panel_canvas, -1, self.figure)
    
    def CleanResults(self):
        self.lc_samples.DeleteAllItems()
    
    def CleanPlot(self):
        self.panel_canvas.DestroyChildren()
        
    def Calculate_dict(self):
        toleran = self.sc_tolerance.GetValue()
        mz_pow = self.sc_mz_power.GetValue()
        int_power = self.sc_intensity.GetValue()
        db = self.cb_database.GetValue()
        if db:
            if not 'Aquasearch_study' in db:
                db = 'Aquasearch_study - ' + db
        dir_ = self.directory_.GetValue()  
        
        if self.folder_sel.IsChecked():
            dic_score = {}
            
            dir_ = dir_.split("\\")[0:-1]
            dir_ = '\\'.join(dir_)
            list_files = os.listdir(dir_)
            
            for file in list_files:
                dir_file = os.path.join(dir_,file)
                
                scores = sm.request_scores(dir_file, tolerance = toleran, 
                                                    shift= 0, mz_power = mz_pow, 
                                                    intensity_power = int_power, dat_b=db)
                
                file = file.split('.')[0]
                dic_score[file] = scores
                
        else:
            dic_score = {}
            
            scores = sm.request_scores(dir_, tolerance = toleran, 
                                                shift= 0, mz_power = mz_pow, 
                                                intensity_power = int_power, dat_b=db)
            
            file = dir_.split("\\")[-1]
            file = file.split('.')[0]
            
            dic_score[file] = scores
        return dic_score
    
    def DownPlot(self, evt):
        plt.savefig('Aquasearch_result.png')

        
if __name__ == "__main__":
            
    MyApp = wx.App()
    aquasearch = AquasearchFrame(None)
    aquasearch.Show()
    MyApp.MainLoop()        