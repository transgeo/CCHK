import pandas as pd
import numpy as np
import sys

sys.setdefaultencoding('utf8')
# Input data
station = pd.read_excel('Interactive_output_HK_ad.xlsx', sheet_name='station')
link = pd.read_excel('Interactive_output_HK_ad.xlsx', sheet_name='line section')
f = open('HKMTRMap.txt', 'w')
# html head body
f.write('<!DOCTYPE html>\r\n'
        '<html lang="en">\r\n'
        '<head>\r\n'
        '   <meta charset="UTF-8">\r\n'
        '   <title>Beijing Metro Map</title>\r\n'
        '</head>\r\n'
        '<body>\r\n'
        '<section class="section" style="padding:0;">\r\n'
        '   <div class="columns is-gapless">\r\n'
        '       <div id="right" class="column" style="border:1px solid black; background-color: white">\r\n'
        '           <div class="Metro Map" id="container" style="width: 100%; height: 100%; ">\r\n'
        '               <svg id="map" xmlns="http://www.w3.org/2000/svg" version="1.1" xlink="http://www.w3.org/1999/xlink"\r\n'
        '                   width="1000" height="1000" xmlns:xlink="http://www.w3.org/1999/xlink"\r\n'
        '                   xmlns:ev="http://www.w3.org/2001/xml-events" style="overflow: hidden; ">\r\n')


# Line section
#  ############ East Rail Line ##########################
# line section
for i in range(link.shape[0]):
    if link['Type'][i] == 2:
        f.write('<path d="M'+bytes(link['Ori_x'][i])+','+bytes(link['Ori_y'][i])+','
                'L'+bytes(link['Des_x'][i])+','+bytes(link['Des_y'][i])+'"\r\n'
                'stroke="'+bytes(link['color'][i])+'" stroke-width="8">\r\n'
                '</path>\r\n')
    else:
        f.write('<path d="M' + bytes(link['Ori_x'][i]) + ',' + bytes(link['Ori_y'][i]) + ','
                'A' + bytes(link['Des_A_r'][i])+','+bytes(link['Des_A_r'][i])+','
                + bytes(link['adt'][i])+','+bytes(int(link['LF'][i]))+','
                + bytes(int(link['SF'][i]))+','
                + bytes(link['Des_x'][i])+','+bytes(link['Des_y'][i])+'"\r\n'
                'fill="none" stroke="'+bytes(link['color'][i])+'" stroke-width="8">\r\n'
                '</path>\r\n')

# station

for i in range(station.shape[0]):
    # text
    f.write('   <text fill="#000"\r\n'
            '       x="'+bytes(station['xx'][i])+'" y="'+bytes(station['yy'][i])+'">\r\n'
            '       <tspan x="'+bytes(station['xx'][i])+'"'
            '               dx="'+bytes(station['label_dx'][i])+'"\r\n'
            '               dy="'+bytes(station['label_dy'][i])+'"\r\n'
            '           style="font-size:12;font-family:SimSun;text-anchor:start;textLength:1;font-weight:normal;"\r\n'
            '           >'+bytes(station['Station'][i])+'</tspan>\r\n'
            '       <tspan x="'+bytes(station['xx'][i])+'"'
            '               dx="'+bytes(station['label_dx'][i])+'"\r\n'
            '               dy="15"\r\n'
            '           style="font-size:10;font-family:Myriad;text-anchor:start;textLength:1;font-weight:normal;"\r\n'
            '           >'+station['Name'][i]+'</tspan>\r\n'
            '   </text>\r\n')
    # Ellipse or transfer station
    if station['Trans'][i] == 0:  # normal station using Ellipse
        f.write('   <ellipse cx="'+bytes(station['xx'][i])+'" cy="'+bytes(station['yy'][i])+'"\r\n'
                '       fill="white"\r\n'
                '       name="'+station['Name'][i]+'" rx="6.5" ry="6.5" stroke="'+bytes(station['color'][i])+'"\r\n'
                '       stroke-width="2.5" sx="100" sy="100">\r\n'
                '   </ellipse>\r\n')
    else:  # transfer station
        f.write('   <svg height="23" version="1.1" viewBox="0 0 105 105" width="23"\r\n'
                '       xmlns="http://www.w3.org/2000/svg"\r\n'
                '       x="'+bytes(station['xx'][i]-10)+'" y="'+bytes(station['yy'][i]-10)+'"\r\n'
                '       name="'+station['Name'][i]+'">\r\n'      
                '       <g id="transfer_station"\r\n>'
                '           <circle cx="52.5" cy="52.5" fill="#ffffff" opacity="1.00" r="46.5" stroke="#231715"\r\n'
                '               stroke-width="12"></circle>\r\n'
                '           <path d=" M 41.04 25.08 C 45.36 23.43 49.98 23.10 54.56 23.33 C 54.46 25.06 54.33 26.78 54.18 28.50 C 48.07 29.49 41.77 31.66 37.79 36.65 C 34.07 42.79 35.78 50.27 34.68 57.03 C 38.42 56.52 42.05 55.42 45.82 55.05 C 40.62 60.86 35.86 67.04 30.95 73.08 C 24.97 67.36 20.01 60.68 14.51 54.51 C 18.24 55.38 22.01 56.08 25.79 56.70 C 25.66 50.84 25.92 44.99 26.63 39.17 C 27.98 32.10 34.23 26.75 41.04 25.08 Z"\r\n'
                '               fill="#231715" opacity="1.00"></path>\r\n'
                '           <path d=" M 72.20 33.89 C 77.70 39.68 82.54 46.05 88.00 51.88 C 84.33 52.04 80.80 51.02 77.23 50.38 C 77.35 55.93 76.91 61.45 76.50 66.98 C 75.98 72.11 72.47 76.44 68.36 79.30 C 62.34 82.73 55.32 84.35 48.40 83.72 C 48.52 82.01 48.64 80.30 48.77 78.59 C 56.66 77.89 64.33 72.81 67.18 65.30 C 68.20 60.37 67.69 55.28 67.82 50.27 C 64.33 50.61 60.97 51.72 57.46 51.93 C 62.48 46.00 67.37 39.97 72.20 33.89 Z"\r\n'                
                '               fill="#231715" opacity="1.00"></path>\r\n'
                '       </g>\r\n'
                '   </svg>')
f.write('</g>')

# html head body
f.write('                       </svg>\r\n'
        '                   </div>\r\n'
        '               </div>\r\n'
        '           </div>\r\n'
        '       </section>\r\n'
        '   </body>\r\n'
        '</html>\r\n')
f.close()
