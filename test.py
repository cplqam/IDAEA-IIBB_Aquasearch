#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import wx
from wx.lib.agw import ultimatelistctrl as ulc

proteins = [('P12345', 'Morsa daurata'),
            ('P23198', 'Pingue beneficius'),
            ('Q43567', 'Roscus sanus')]


class MyFrame(wx.Frame):
    def __init__(self, *args, **kwds):
        wx.Frame.__init__(self, *args, **kwds)
        self.SetSize((400, 300))
        self.SetTitle("test_frame")
        
        self.panel = wx.Panel(self, wx.ID_ANY)
        sizer = wx.BoxSizer(wx.VERTICAL)
        
        self.ulist = ulc.UltimateListCtrl(self.panel, wx.ID_ANY, agwStyle=wx.LC_HRULES | wx.LC_REPORT | wx.LC_VRULES)
        self.ulist.InsertColumn(col=0, heading="CODE", format=wx.LIST_FORMAT_LEFT, width=-1)
        self.ulist.InsertColumn(col=1, heading= "NAME", format=wx.LIST_FORMAT_LEFT, width=200)
        sizer.Add(self.ulist, 1, wx.EXPAND, 0)
        
        self.panel.SetSizer(sizer)
        self.Layout()

        font = wx.Font(wx.FontInfo(9).Italic())

        for idx, (code, name) in enumerate(proteins):
            self.ulist.InsertStringItem(idx, code, 0)
            self.ulist.SetStringItem(idx, 1, name)
            item = self.ulist.GetItem(idx, 1)
            item.SetMask(ulc.ULC_MASK_FONT)
            item.SetFont(font)
            self.ulist.SetItem(item)

if __name__ == "__main__":
    app = wx.App()
    frame = MyFrame(None, wx.ID_ANY, "")
    frame.Show()
    app.MainLoop()
