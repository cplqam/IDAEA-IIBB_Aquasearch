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
from main_aquasearch import MainAquasearch
from results_aqua import ResultsAqua
from results_peptides import ResultsPeptides
from wx.lib.agw import ultimatelistctrl as ulc
import score_matchms as sm
import load_archives
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
class AquasearchFrame(MainAquasearch):

    def __init__(self, *args, **kwds):

        MainAquasearch.__init__(self, *args, **kwds)

        self.dic_s = {}
        self.Bind(wx.EVT_BUTTON, self.on_bt_query, self.button_run_query)
        self.Bind(wx.EVT_BUTTON, self.on_bt_plot, self.button_pca)
        self.Bind(wx.EVT_BUTTON, self.on_bt_1pt_mix, self.btn_1lbl_mix_run)
        self.Bind(wx.EVT_BUTTON, self.on_bt_1pt_std, self.btn_1lbl_std_run)
        self.Bind(wx.EVT_BUTTON, self.on_bt_2pt_mix, self.btn_2lbl_mix_run)
        self.Bind(wx.EVT_BUTTON, self.on_bt_2pt_std, self.btn_2lbl_std_run)
        self.Bind(wx.EVT_BUTTON, self.down_plot, self.button_download)
        self.Bind(wx.EVT_BUTTON, self.on_bt_del_prot, self.btn_del_prot)
        self.Bind(wx.EVT_BUTTON, self.on_bt_del_sam, self.btn_del_sam)

        self.Bind(wx.EVT_COMBOBOX, self.on_selection_protein, self.cb_database_delete_prot)
        self.Bind(wx.EVT_COMBOBOX, self.on_selection_sample, self.cb_database_delete_sam)
        self.Bind(wx.EVT_COMBOBOX, self.new_db, self.cb_database_ne)

        self.Bind(wx.EVT_SPINCTRLDOUBLE, self.components_sel, self.sc_number_PCs)

        self.resized = True
        self.Bind(wx.EVT_SIZE,self.on_size)
        self.Bind(wx.EVT_IDLE,self.on_idle)

        dbs_ =sr.db_request()
        self.cb_database.Clear()
        for cb_option in dbs_:
            self.cb_database.Append(cb_option)

        self.cb_database_delete_prot.Clear()
        for cb_option in dbs_:
            self.cb_database_delete_prot.Append(cb_option)
            
        self.cb_database_delete_sam.Clear()
        for cb_option in dbs_:
            self.cb_database_delete_sam.Append(cb_option)
        
        dbs_.append('New database (select to create a new database)')
        self.cb_database_ne.Clear()
        for cb_option in dbs_:
            self.cb_database_ne.Append(cb_option)
        
        self.text_tolerance.SetToolTip(hs.tolerance)
        self.text_db_query.SetToolTip(hs.database)
        self.text_adv_det.SetToolTip(hs.advanced_details)
        self.txt_db_ne.SetToolTip(hs.database_ne)
        self.text_prot_results.SetToolTip(hs.prot_results)
        self.text_pcs.SetToolTip(hs.pcs)
        self.text_pcx.SetToolTip(hs.pcx)
        self.txt_pcy.SetToolTip(hs.pcy)
        self.txt_pcz.SetToolTip(hs.pcz)
        self.text_title.SetToolTip(hs.title)
        self.txt_size.SetToolTip(hs.size)
        
        self.txt_maldi_ind_mix.SetToolTip(hs.maldi_ind_mix)
        self.txt_pep_ind_mix.SetToolTip(hs.pep_ind_mix)
        self.txt_prot_ind_mix.SetToolTip(hs.prot_ind_mix)
        self.txt_ppm_ind_mix.SetToolTip(hs.ppm_ind_mix)
        self.txt_uniq_ind_mix.SetToolTip(hs.unique_ind_mix)
        self.txt_code_ind_mix.SetToolTip(hs.code_ind_mix)
        self.txt_name_ind_mix.SetToolTip(hs.name_ind_mix)
        self.txt_rank_ind_mix.SetToolTip(hs.rank_ind_mix)
        self.txt_maldi_indiv_std.SetToolTip(hs.maldi_indiv_std)
        self.label_uniq_indiv_std.SetToolTip(hs.unique_ind_std)
        self.label_ppm_indiv_std.SetToolTip(hs.ppm_ind_std)
        self.txt_code_indiv_std.SetToolTip(hs.code_ind_std)
        self.txt_name__indiv_std.SetToolTip(hs.sample_ind_std)
        self.txt_maldi_gru_mix.SetToolTip(hs.maldi_ind_mix)
        self.txt_peptide_gru_mix.SetToolTip(hs.pep_ind_mix)
        self.txt_prot_gru_mix.SetToolTip(hs.prot_ind_mix)
        self.txt_ppm_gru_mix.SetToolTip(hs.ppm_ind_mix)
        self.txt_uniq_gru_mix.SetToolTip(hs.unique_ind_mix)
        self.txt_code1_gru_mix.SetToolTip(hs.code1_gru_mix)
        self.txt_name_gru_mix.SetToolTip(hs.name_gru_mix)
        self.txt_sam_gru_mix.SetToolTip(hs.name_ind_mix)
        self.txt_rank_gru_mix.SetToolTip(hs.rank_ind_mix)
        self.txt_maldi_gru_std.SetToolTip(hs.maldi_indiv_std)
        self.txt_uniq_gru_std.SetToolTip(hs.unique_ind_std)
        self.txt_ppm_gru_std.SetToolTip(hs.ppm_ind_std)
        self.txt_code_gru_std.SetToolTip(hs.code_gru_std)
        self.txt_name_gru_std.SetToolTip(hs.name_gru_std)
        self.txt_sam_gru_std.SetToolTip(hs.name_ind_mix)
        
        self.txt_del_db1.SetToolTip(hs.del_db1)
        self.txt_del_db2.SetToolTip(hs.del_db2)
        self.txt_del_prot.SetToolTip(hs.del_prot)
        self.txt_del_samp.SetToolTip(hs.del_sam)
        
    def new_db(self, evt):
        try: 
            dbs_ = sr.db_request()
            dbs_.append('New database (select to create a new database)')
            db = self.cb_database_ne.GetSelection()
            db = dbs_[db]
            if db == "New database (select to create a new database)":
                new_db = ''
                while not new_db:
                    dlg = wx.TextEntryDialog(self,'New database', 'Database name', "",
                                             wx.OK|wx.CANCEL)
                    answ = dlg.ShowModal()
                    if answ == wx.ID_OK:
                        new_db = dlg.GetValue()
                        if not new_db:
                            wx.MessageBox('Error: Introduce the name of the new database',
                                          'Info', wx.OK | wx.ICON_INFORMATION)
                    
                    elif answ == wx.ID_CANCEL:  
                        break
                        
                if not new_db:
                    dlg.Close()
                else:
                    sr.create_db_by_user(new_db)
                    new_list = sr.db_request()
                        
                    self.cb_database_delete_prot.Clear()
                    for cb_option in new_list:
                        self.cb_database_delete_prot.Append(cb_option)
                        
                    self.cb_database_delete_sam.Clear()
                    for cb_option in new_list:
                        self.cb_database_delete_sam.Append(cb_option)
                
                    self.cb_database.Clear()
                    for cb_option in new_list:
                        self.cb_database.Append(cb_option)
                    
                    new_list.append("New database (select to create a new database)")
                    self.cb_database_ne.Clear()
                    for cb_option in new_list:
                        self.cb_database_ne.Append(cb_option)
                
        except sqlite3.OperationalError: 
            wx.MessageBox('Error: No database selected', 'Info', wx.OK | wx.ICON_INFORMATION)
                
    
    def on_bt_del_prot(self, evt):
        try:
            db_ = self.cb_database_delete_prot.GetValue()
            protein = self.cb_delete_prot.GetValue()
            if db_ and protein:
            
                dlg = wx.MessageDialog(None, 'Are you sure you want to delete the protein "' + protein + '"?',
                                       'Confirm Delete', wx.YES_NO | wx.ICON_QUESTION)
                dlg.SetYesNoLabels('Yes', 'No')
                answ = dlg.ShowModal()
                if answ == wx.ID_YES:
                    wx.BeginBusyCursor(wx.Cursor(wx.CURSOR_ARROWWAIT))
                    dbe.protein_elimination(db_, protein)
                    self.on_selection_protein(self.cb_database_delete_prot)
                    wx.EndBusyCursor()
                    wx.MessageBox('Process has been successfully completed !!', 'Info',
                                  wx.OK | wx.ICON_INFORMATION)
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
                    self.on_selection_sample(self.cb_database_delete_sam)
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
    
    def on_selection_protein(self, evt):
        dbs_option = sr.db_request()
        print(dbs_option)
        db_chosen = self.cb_database_delete_prot.GetSelection()
        db_chosen = dbs_option[db_chosen]
        try:
            if 'Aquasearch_study' not in db_chosen:
                db_chosen = 'Aquasearch_study - ' + db_chosen
        
            print(db_chosen)
            table_ = pd.DataFrame(sr.table_download(db_chosen, 'Protein_codes'))
           
            self.cb_delete_prot.Clear()
            for cod in table_.iloc[:,1]:
                self.cb_delete_prot.Append(cod)
            
        except sqlite3.OperationalError:
            pass
    
    def on_selection_sample(self, evt):
        dbs_option =sr.db_request()
        db_chosen = self.cb_database_delete_sam.GetSelection()
        db_chosen = dbs_option[db_chosen]
        try:
            if 'Aquasearch_study' not in db_chosen:
                db_chosen = 'Aquasearch_study - ' + db_chosen
            
            table_ = pd.DataFrame(sr.table_download(db_chosen, 'Spectra_table'))
            table_ = table_.iloc[:,4].unique()
            table_.sort()
        
            self.cb_delete_sam.Clear()
            for cod in table_:
                self.cb_delete_sam.Append(cod)
        except sqlite3.OperationalError: 
            pass
    
    def components_sel(self, evt):
        n_pc = int(self.sc_number_PCs.GetValue())
        
        self.PC_x.Clear()
        self.PC_y.Clear()
        self.PC_z.Clear()
        
        self.PC_z.Append('None')
        for num in range(n_pc):
            self.PC_x.Append(str(num+1))
            self.PC_y.Append(str(num+1))
            self.PC_z.Append(str(num+1))
    
    def on_size(self, evt):
        self.resized = True
        evt.Skip()
        
    def on_idle(self, evt):
        if self.resized: 
            self.clean_plot()
            self.figure.set_size_inches(self.panel_canvas.Size[0]/105, self.panel_canvas.Size[1]/105)
            self.canvas = FigureCanvas(self.panel_canvas, -1,  self.figure)
            
            self.resized = False
    
    def on_bt_query(self, evt):
        try:
            wx.BeginBusyCursor(wx.Cursor(wx.CURSOR_ARROWWAIT))
            self.clean_results()
        
            self.dic_s = self.calculate_dict()
        
            for numero in self.dic_s.keys():
                self.lc_samples.InsertItem(1,numero)
            
            self.Bind(wx.EVT_LIST_ITEM_ACTIVATED, lambda evt, temp=self.dic_s: self.sample_display(evt, temp))
            wx.EndBusyCursor()
            
        except FileNotFoundError:
            wx.MessageBox('Error: No file pathway selected', 'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
        except sqlite3.OperationalError: 
            wx.MessageBox('Error: No database selected', 'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
        except UnicodeDecodeError:
            wx.MessageBox('Error: A txt file with MALDI-TOF results have to be selected',
                          'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
        
    def on_bt_plot(self, evt):
        try: 
            self.clean_plot()
            
            wx.BeginBusyCursor(wx.Cursor(wx.CURSOR_ARROWWAIT))
            if not self.dic_s:
                self.dic_s = self.calculate_dict()
                wx.EndBusyCursor()
        
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
            wx.EndBusyCursor()
                
        except sqlite3.OperationalError: 
            wx.MessageBox('Error: No database selected', 'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
        except FileNotFoundError:
            wx.MessageBox('Error: No file pathway selected', 'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
        except ValueError as err:
            if "invalid literal for int() with base 10: ''" in str(err):
                wx.MessageBox('Error: Select all PC correctly', 'Info', wx.OK | wx.ICON_INFORMATION)
                wx.EndBusyCursor()
            elif "Input X contains NaN" in str(err):
                wx.MessageBox("""Error: PCA analysis can not be run if any sample has no identifications. 
                              Remove this/these sample/s and repeat the analysis""",
                              'Info', wx.OK | wx.ICON_INFORMATION)
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
            
            rank = self.sc_1prot_rank.GetValue()
            sample_name = self.txt_1prot_mix_name.GetValue()
            code = self.txt_1prot_mix_code.GetValue()
            
            if self.all_prot_ind_mix.IsChecked():
                if not sample_name:
                    raise ValueError
                    
                db_data = load_archives.db_table_request(db)
                db_id = list(db_data.keys())                    
                test_pdmm = pdmm.txt_complete(maldi, peptide, protein, rank, ppm, uni_pep, db)
                
                for prot_id in db_id:
                    print(prot_id)
                    protein_code = sr.consulta_prot_id('Aquasearch_study', prot_id)
                    if not '(' in protein_code:
                        ps.fill_table(protein_code, test_pdmm, sample_name, db)
                        print(protein_code)
                    else:
                        pass
                
            else:
                if not code or not sample_name:
                    raise ValueError
        
                test_pdmm = pdmm.txt_complete(maldi, peptide, protein, rank, ppm, uni_pep, db)

                ps.fill_table(code, test_pdmm, sample_name, db)
            
            wx.EndBusyCursor()
            wx.MessageBox('Process has been successfully completed !!', 'Info', wx.OK | wx.ICON_INFORMATION)
        except sqlite3.OperationalError: 
            wx.MessageBox('Error: No database selected or the selected database is empty and "all inividual protein codes" is selected',
                          'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
        except FileNotFoundError:
            wx.MessageBox('Error: No file pathway selected', 'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
        except ValueError:
            wx.MessageBox('Error: Select correctly all the parameters', 'Info', wx.OK | wx.ICON_INFORMATION)
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
            wx.MessageBox('Error: Select correctly all the parameters', 'Info', wx.OK | wx.ICON_INFORMATION)
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
            
            if self.all_proteins_gr.IsChecked():
                if not sample_name:
                    raise ValueError
                    
                db_data = load_archives.db_table_request(db)
                db_id = list(db_data.keys())                    
                test_pdmm = pdmm.txt_complete(maldi, peptide, protein, rank, ppm, uni_pep, db)
                
                for prot_id in db_id:
                    protein_code = sr.consulta_prot_id('Aquasearch_study', prot_id)
                    if '(' in protein_code:
                        print(protein_code)
                        protein_code = protein_code.split('_(')
                        new_entry = protein_code[0]
                        print(new_entry)
                        protein_codes_joined = protein_code[1]
                        protein_codes_joined = protein_codes_joined.replace(')', '').split(';')
                        print(protein_codes_joined)
                        
                        ps.table_union(new_entry, protein_codes_joined[0], protein_codes_joined[1], test_pdmm, sample_name, db)
                    else:
                        pass
            else:
            
                if not code1 or not code2 or not new_entry or not sample_name:
                    raise ValueError
        
                test_pdmm = pdmm.txt_complete(maldi, peptide, protein, rank, ppm, uni_pep, db)
                ps.table_union(new_entry, code1,code2, test_pdmm, sample_name, db)
            
            wx.EndBusyCursor()
            wx.MessageBox('Process has been successfully completed !!', 'Info', wx.OK | wx.ICON_INFORMATION)
        except sqlite3.OperationalError: 
            wx.MessageBox('Error: No database selected or the selected database is empty and "all grouped protein codes" is selected',
                          'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
        except FileNotFoundError:
            wx.MessageBox('Error: No file pathway selected', 'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
        except ValueError:
            wx.MessageBox('Error: Select correctly all the parameters', 'Info', wx.OK | wx.ICON_INFORMATION)
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
            wx.MessageBox('Error: Select correctly all the parameters', 'Info', wx.OK | wx.ICON_INFORMATION)
            wx.EndBusyCursor()
            
    def draw_3d(self, evt, data_matrix, var, x_value, y_value, z_value, title, size_point, sample_names):
    
        self.figure = plt.figure(figsize=(self.panel_canvas.Size[0]/105, self.panel_canvas.Size[1]/105))
        ax_3d = self.figure.add_subplot(projection='3d')
        ax_3d.scatter(data_matrix[:,x_value-1], data_matrix[:,y_value-1], data_matrix[:,z_value-1], s = size_point)
        
        for i, txt in enumerate(sample_names): 
            ax_3d.text3D(data_matrix[i,x_value-1]+0.015, data_matrix[i,y_value-1]+0.015, data_matrix[i,z_value-1]+0.015,txt)
            
        ax_3d.set_xlabel('PC' + str(x_value) + ' (' + str(round(var[x_value-1]*100, 2)) + '%)')
        ax_3d.set_ylabel('PC' + str(y_value) + ' (' + str(round(var[y_value-1]*100, 2)) + '%)')
        ax_3d.set_zlabel('PC' + str(z_value) + ' (' + str(round(var[z_value-1]*100, 2)) + '%)')
        plt.title(title)
        self.canvas = FigureCanvas(self.panel_canvas, -1, self.figure)
    
    def draw_2d(self, evt, data_matrix, var, x_value, y_value, title, size_point, sample_names):
        self.figure = plt.figure(figsize=(self.panel_canvas.Size[0]/105, self.panel_canvas.Size[1]/105))
        ax_2d = self.figure.add_subplot()
        ax_2d.scatter(data_matrix[:,x_value-1], data_matrix[:,y_value-1], s = size_point)
        
        for i, txt in enumerate(sample_names):
            ax_2d.annotate(txt, (data_matrix[i,x_value-1]+0.01, data_matrix[i,y_value-1]+0.01))
            
        ax_2d.set_xlabel('PC' + str(x_value) + ' (' + str(round(var[x_value-1]*100, 2)) + '%)')
        ax_2d.set_ylabel('PC' + str(y_value) + ' (' + str(round(var[y_value-1]*100, 2)) + '%)')
        plt.title(title)
        self.canvas = FigureCanvas(self.panel_canvas, -1, self.figure)
    
    def clean_results(self):
        self.lc_samples.DeleteAllItems()
    
    def clean_plot(self):
        self.panel_canvas.DestroyChildren()
        
    def calculate_dict(self):
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
    
    def down_plot(self, evt):
        with wx.FileDialog(self, "Save PNG file", wildcard="PNG files (*.png)|*.png",
                       style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as fileDialog:

            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return     # the user changed their mind

            # save the current contents in the file
            pathname = fileDialog.GetPath()
            try:
                plt.savefig(pathname)
            except IOError:
                wx.LogError("Cannot save current data in file '%s'." % pathname)
    
    def fbbCallback(self, evt):
        pass

    def sample_display(self, event, dic_s):
        
        name = event.GetText()
        score = dic_s[name]
        
        self.results_sample = ResultsAqua(name, None, wx.ID_ANY, "")
        
        new_proteins = []
        font = wx.Font(8, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_ITALIC, wx.FONTWEIGHT_NORMAL)
        for proteins in zip(score.keys(), score.values()):
            tup = (proteins[0], proteins[1][0][0], proteins[1][0][1], str(proteins[1][0][2]), str(proteins[1][0][3]), 
                    str(proteins[1][0][4]))
            new_proteins.append(tup)
        
        for idx, (code, name_p, name_o, sco_, pep_t, pep_u) in enumerate(new_proteins):
            self.results_sample.ulist.InsertStringItem(idx, code, 0)
            self.results_sample.ulist.SetStringItem(idx, 1, name_p)
            self.results_sample.ulist.SetStringItem(idx, 2, name_o)
            self.results_sample.ulist.SetStringItem(idx, 3, sco_)
            self.results_sample.ulist.SetStringItem(idx, 4, pep_t)
            self.results_sample.ulist.SetStringItem(idx, 5, pep_u)
            item = self.results_sample.ulist.GetItem(idx, 2)
            item.SetMask(ulc.ULC_MASK_FONT)
            item.SetFont(font)
            self.results_sample.ulist.SetItem(item)
        
        self.results_sample.Bind(wx.EVT_LIST_ITEM_ACTIVATED, lambda evt, temp=score: self.peptides_display(evt, temp))
        
        self.results_sample.Show()
    
    
    def peptides_display(self, event, protein_inf):
        protein = event.GetText()
        pep_information = protein_inf[protein][1]
        
        self.results_protein = ResultsPeptides(protein, None, wx.ID_ANY, "")
        
        for pos in range(pep_information.shape[0]):
            index = self.results_protein.list_results.InsertItem(self.results_protein.list_results.GetItemCount(), 
                                                                 pep_information.iloc[pos,0])
            self.results_protein.list_results.SetItem(index, 1, str(pep_information.iloc[pos,1]))
            self.results_protein.list_results.SetItem(index, 2, str(pep_information.iloc[pos,2]))
        
        self.results_protein.Show()
        
        
    def onclick(self, evt):
        # EXITS APPLICATION ON CLICKING EXIT BUTTON
        self.Close()
        
if __name__ == "__main__":
            
    MyApp = wx.App()
    aquasearch = AquasearchFrame(None)
    aquasearch.Show()
    MyApp.MainLoop()        