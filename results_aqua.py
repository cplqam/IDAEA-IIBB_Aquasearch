#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
#
# generated by wxGlade 1.0.4 on Wed Feb 22 13:08:05 2023
#

import wx
from wx.lib.agw import ultimatelistctrl as ulc
# begin wxGlade: dependencies
# end wxGlade

# begin wxGlade: extracode
# end wxGlade


class ResultsAqua(wx.Frame):
    def __init__(self, name, *args, **kwds):
        # begin wxGlade: results_aqua.__init__
        kwds["style"] = kwds.get("style", 0) | wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwds)
        self.SetSize((700, 600))
        self.SetTitle("Aquasearch results")
        _icon = wx.NullIcon
        _icon.CopyFromBitmap(wx.Bitmap("Icon.png", wx.BITMAP_TYPE_ANY))
        self.SetIcon(_icon)

        self.panel_1 = wx.Panel(self, wx.ID_ANY)

        sizer_1 = wx.BoxSizer(wx.VERTICAL)

        self.panel_2 = wx.Panel(self.panel_1, wx.ID_ANY, style=wx.BORDER_STATIC)
        self.panel_2.SetMinSize((-1, 40))
        self.panel_2.SetBackgroundColour(wx.Colour(255, 255, 0))
        sizer_1.Add(self.panel_2, 0, wx.ALL | wx.EXPAND, 5)

        sizer_2 = wx.BoxSizer(wx.HORIZONTAL)

        self.panel_4 = wx.Panel(self.panel_2, wx.ID_ANY)
        sizer_2.Add(self.panel_4, 1, wx.ALIGN_CENTER_VERTICAL, 0)

        sizer_3 = wx.BoxSizer(wx.VERTICAL)

        self.label_1 = wx.StaticText(self.panel_4, wx.ID_ANY, "Results of sample " + name)
        self.label_1.SetFont(wx.Font(14, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_ITALIC, wx.FONTWEIGHT_NORMAL, 0, "Segoe UI"))
        sizer_3.Add(self.label_1, 0, wx.ALIGN_CENTER_HORIZONTAL, 0)

        self.panel_3 = wx.Panel(self.panel_1, wx.ID_ANY, style=wx.BORDER_STATIC)
        self.panel_3.SetBackgroundColour(wx.Colour(255, 255, 255))
        sizer_1.Add(self.panel_3, 1, wx.EXPAND | wx.LEFT | wx.RIGHT, 5)

        sizer_5 = wx.BoxSizer(wx.HORIZONTAL)

        self.ulist = ulc.UltimateListCtrl(self.panel_3, wx.ID_ANY, agwStyle=wx.LC_HRULES | wx.LC_REPORT | wx.LC_VRULES)
        self.ulist.InsertColumn(col=0, heading="Protein code", format=wx.LIST_FORMAT_LEFT, width=110)
        self.ulist.InsertColumn(col=1, heading="Protein name", format=wx.LIST_FORMAT_LEFT, width=130)
        self.ulist.InsertColumn(col=2, heading="Organism", format=wx.LIST_FORMAT_LEFT, width=150)
        self.ulist.InsertColumn(col=3, heading="Score", format=wx.LIST_FORMAT_LEFT, width=80)
        self.ulist.InsertColumn(col=4, heading="Nº of peptides", format=wx.LIST_FORMAT_LEFT, width=90)
        self.ulist.InsertColumn(col=5, heading="Unique peptides", format=wx.LIST_FORMAT_LEFT, width=100)
        
        sizer_5.Add(self.ulist, 1, wx.EXPAND, 0)

        self.panel_5 = wx.Panel(self.panel_1, wx.ID_ANY)
        self.panel_5.SetMinSize((-1, 30))
        sizer_1.Add(self.panel_5, 0, wx.EXPAND, 0)

        sizer_4 = wx.BoxSizer(wx.HORIZONTAL)

        self.btn_export = wx.Button(self.panel_5, wx.ID_ANY, "Export")
        sizer_4.Add(self.btn_export, 0, wx.ALL, 5)

        self.panel_7 = wx.Panel(self.panel_5, wx.ID_ANY)
        sizer_4.Add(self.panel_7, 1, wx.EXPAND, 0)

        self.btn_close = wx.Button(self.panel_5, wx.ID_ANY, "Close")
        self.Bind(wx.EVT_BUTTON, self.onclick, self.btn_close)
        sizer_4.Add(self.btn_close, 0, wx.ALL, 5)

        self.panel_5.SetSizer(sizer_4)

        self.panel_3.SetSizer(sizer_5)

        self.panel_4.SetSizer(sizer_3)

        self.panel_2.SetSizer(sizer_2)

        self.panel_1.SetSizer(sizer_1)

        self.Layout()
        # end wxGlade
    def onclick(self, e):
        # EXITS APPLICATION ON CLICKING EXIT BUTTON
        self.Close()
    

# end of class results_aqua

class MyApp(wx.App):
    def OnInit(self):
        self.Results = ResultsAqua(self, wx.ID_ANY, "")
        self.SetTopWindow(self.Results)
        self.Results.Show()
        return True

# end of class MyApp

if __name__ == "__main__":
    app = MyApp(0)
    app.MainLoop()
