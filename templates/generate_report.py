#!/usr/bin/env python3

import sys
import re
import numpy as np
import pandas as pd
import base64
from math import modf
import datetime
pd.set_option('max_colwidth',200)

if len(sys.argv) == 1:
    sys.argv=["generate_report.py","$base","$orig","$dup_marks","$snpmisspng","$indmisspng","$maf_png","$maf_txt","$hwe_png","$hwe_qq_png","$qc1_out","$badsex","$pca_png","$fail_IBD","$miss_het_png","$fail_het","$re_maf","$re_hwe","$newbim","$newfam","$qc2_out"]

sample_name = sys.argv[1]
basic = sys.argv[2]
dups = sys.argv[3]

snpmisspng = sys.argv[4]
indmisspng = sys.argv[5]
maf_png = sys.argv[6]
maf_txt = sys.argv[7]
hwe_png = sys.argv[8]
hwe_qq_png = sys.argv[9]
qc1_out = sys.argv[10]

badsex = sys.argv[11]
pca_png = sys.argv[12]
fail_IBD = sys.argv[13]
miss_het_png = sys.argv[14]
fail_het = sys.argv[15]
qc2_out = sys.argv[20]

re_maf = sys.argv[16]
re_hwe = sys.argv[17]

newbim = sys.argv[18]
newfam = sys.argv[19]

today = datetime.date.today()

GEN_HTML = sample_name+'-GWAS-QC_report.html'

def dftohtml(file):
    if type(file) == str:
        df = pd.read_table(file)
    else:
        df = file
    a = df.to_html(index=False)
    b = a.split('\\n')
    if b[3] == '      <th>Unnamed: 0</th>':
        b[3] ='       <th></th>'
    c = ''.join(b[1:-1])
    return c

def floattoint(a):
    if modf(a)[0] == 0:
        a = str(int(a))
    else:
        a = str(round(a, 6))
    return a

def pngtobase64(png):
    f = open(png,"rb")
    base64_data = base64.b64encode(f.read()).decode("utf-8")
    return base64_data

def table_html(df,module_id,name):
    html = """
<div class="module">
    <h2 id=%s>
        %s
    </h2>
    <table class="dataintable">%s</table>
</div>
"""%(module_id,name,df)
    return html

def png_html(png_code,module_id,name):
    html = """
<div class="module">
    <h2 id=%s>
        %s
    </h2>
    <p>
        <img class="indented" src="data:image/png;base64,%s" alt="01pie" width="600" height="600"/>
    </p>
</div>
"""%(module_id,name,png_code)
    return html

basic_df = pd.read_table(basic,sep=' ',names=['Num of Input Markers','Input Files'])
basic_df['Num of Output Markers'] = [len(open(newbim).readlines()),len(open(newfam).readlines())]
basic_df['Output Files'] = [newbim,newfam]
basic_df['Removed'] = [basic_df['Num of Input Markers'][0] - basic_df['Num of Output Markers'][0],basic_df['Num of Input Markers'][1] - basic_df['Num of Output Markers'][1]]

dup_df = pd.DataFrame(columns=['Num of Dupmarkers','File of Dupmarkers'])
dup_df.loc[0] = [len(open(dups).readlines()),dups]
dup_html = dftohtml(dup_df)

maf_df = pd.read_table(maf_txt)
maf_df.columns = ['MAF bins','Num of Markers']

i = 0
with open(badsex,'r') as f:
    for line in f:
        line_sp = line.split()
        if i == 0:
            badsex_df = pd.DataFrame(columns=line_sp)
            i += 1
        else:
            badsex_df.loc[i] = line_sp
            i += 1

navigation = """
    <body>
        <div class="header">
            <div id="header_title">
                GWAS QC Report
            </div>
            <div id="header_filename"> """+str(today)+""" <br/> """+sample_name+"""
            </div>
        </div>        
        <div class="summary">
            <h2>Contents</h2>
            <ul>
                <li>
                    <a href="#M0">0.   Input & Output</a>
                </li>
                <li>
                    <a href="#M1">1.   Remove duplicated Markers</a>
                </li>
                <li>
                    <a href="#M2">2.1  Markers GT Missing</a>
                </li>
                <li>
                    <a href="#M3">2.2  Individuals GT Missing</a>
                </li>
                <li>
                    <a href="#M4">2.3  MAF</a>
                </li>
                <li>
                    <a href="#M5">2.4  HWE</a>
                </li>
                <li>
                    <a href="#M6">2.5  Phase 2 summary</a>
                </li>
                <li>
                    <a href="#M7">3.1  Individuals with dadsex</a>
                </li>
                <li>
                    <a href="#M8">3.2  PCA</a>
                </li>
                <li>
                    <a href="#M9">3.3  Individual fail IBD</a>
                </li>
                <li>
                    <a href="#M10">3.4  Check het</a>
                </li>
                <li>
                    <a href="#M11">3.5 Phase 3 summary</a>
                </li>
                <li>
                    <a href="#M12">4.1  Review MAF</a>
                </li>
                <li>
                    <a href="#M13">4.2  Review HWE</a>
                </li>
            </ul>
        </div>
        <div class="main">
        """

tail =  """</div>
        <div class="footer">
            Produced by Jianhua Wang
        </div>
    </body>
</html>"""

header = """<html>
    <head>
        <meta charset="utf-8"> 
        <title> GWAS-QC Report</title>
        <style type="text/css">
            @media screen {
            div.summary {
            width: 18em;
            position:fixed;
            top: 3em;
            margin:1em 0 0 1em;
            }
            
            div.main {
            display:block;
            position:absolute;
            overflow:auto;
            height:auto;
            width:auto;
            top:4.5em;
            bottom:2.3em;
            left:18em;
            right:0;
            border-left: 1px solid #CCC;
            padding:0 0 0 1em;
            background-color: white;
            z-index:1;
            }
            
            div.header {
            background-color: #EEE;
            border:0;
            margin:0;
            padding: 0.5em;
            font-size: 200%;
            font-weight: bold;
            position:fixed;
            width:100%;
            top:0;
            left:0;
            z-index:2;
            }
        
            div.footer {
            background-color: #EEE;
            border:0;
            margin:0;
            padding:0.5em;
            height: 1.3em;
            overflow:hidden;
            font-size: 100%;
            font-weight: bold;
            position:fixed;
            bottom:0;
            width:100%;
            z-index:2;
            }
            
            img.indented {
            margin-left: 3em;
            }
            }
            
            @media print {
            img {
                max-width:100% !important;
                page-break-inside: avoid;
            }
            h2, h3 {
                page-break-after: avoid;
            }
            div.header {
                background-color: #FFF;
            }
            
            }
            
            body {    
            font-family: sans-serif;   
            color: #000;   
            background-color: #FFF;
            border: 0;
            margin: 0;
            padding: 0;
            }
            
            div.header {
            border:0;
            margin:0;
            padding: 0.5em;
            font-size: 200%;
            font-weight: bold;
            width:100%;
            }    
            
            #header_title {
            display:inline-block;
            float:left;
            clear:left;
            }
            #header_filename {
            display:inline-block;
            float:right;
            clear:right;
            font-size: 50%;
            margin-right:2em;
            text-align: right;
            }
        
            div.header h3 {
            font-size: 50%;
            margin-bottom: 0;
            }
                
            div.main {
            background-color: white;
            }
                
            div.module {
            padding-bottom:1.5em;
            padding-top:1.5em;
            }
                
            div.footer {
            background-color: #EEE;
            border:0;
            margin:0;
            padding: 0.5em;
            font-size: 100%;
            font-weight: bold;
            width:100%;
            }
        
            div.summary ul {
            padding-left:0;
            list-style-type:none;
            }

            ul a{
            display: block;
            width: 300px;
            height: 40px;
            line-height: 40px;
            background: #f8f8f8;
            margin-bottom: -1px;
            }

            a:link, a:visited {text-decoration: none;}
            a:hover, a:active {text-decoration: underline;}

            ul li:last-child a{
            margin-bottom: 0;
            }

            ul .active a, ul a:hover{
            background: #e7e7e7;
            }
                
            h2 {
            color: #800000;
            padding-bottom: 0;
            margin-bottom: 0;
            clear:left;
            }
        
            table.dataintable {
            margin-top:15px;
            border-collapse:collapse;
            border:1px solid #aaa;
            width:50%;
            }
            table.dataintable th {
            vertical-align:baseline;
            padding:5px 15px 5px 6px;
            background-color:#3F3F3F;
            border:1px solid #3F3F3F;
            text-align:left;
            color:#fff;
            }
            table.dataintable td {
            vertical-align:text-top;
            padding:6px 15px 6px 6px;
            border:1px solid #aaa;
            }
            table.dataintable tr:nth-child(odd) {
            background-color:#F5F5F5;
            }
            table.dataintable tr:nth-child(even) {
            background-color:#fff;
            }
        
            img {
            padding-top: 0;
            margin-top: 0;
            border-top: 0;
            }
        
            
            p {
            padding-top: 0;
            margin-top: 0;
            }
        </style>
    </head>
    """

message = header+navigation+\
table_html(dftohtml(basic_df),'"M0"','0. Input & Output')+\
table_html(dftohtml(dup_df),'"M1"','1. Remove duplicated Markers')+\
png_html(pngtobase64(snpmisspng),'"M2"','2.1 Markers GT Missing')+\
png_html(pngtobase64(indmisspng),'"M3"','2.2 Individuals GT Missing')+\
png_html(pngtobase64(maf_png),'"M4"','2.3 MAF')+\
table_html(dftohtml(maf_df),'"M4.1"','MAF counts')+\
png_html(pngtobase64(hwe_png),'"M5"','2.4 HWE')+\
png_html(pngtobase64(hwe_qq_png),'"M5.1"','HWE QQ plot')+\
table_html(dftohtml(qc1_out),'"M6"','2.5 Phase 2 summary')+\
table_html(dftohtml(badsex_df),'"M7"','3.1 Individuals with dadsex')+\
png_html(pngtobase64(pca_png),'"M8"','3.2 PCA')+\
table_html(dftohtml(fail_IBD),'"M9"','3.3 Individual fail IBD')+\
png_html(pngtobase64(miss_het_png),'"M10"','3.4 Check het')+\
table_html(dftohtml(fail_het),'"M10.1"','Individual fail HET')+\
table_html(dftohtml(qc2_out),'"M11"','3.5 Phase 3 summary')+\
png_html(pngtobase64(re_maf),'"M12"','4.1 Review MAF')+\
png_html(pngtobase64(re_hwe),'"M13"','4.2 Review HWE')+\
tail

f = open(GEN_HTML,'w')
f.write(message)
f.close()