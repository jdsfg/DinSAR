import sys, math
import numpy as np

def fix_prm(prm_file, xml_file, led_file):
    import xml.etree.ElementTree as ET
    tree = ET.parse(xml_file)
    root = tree.getroot()
    lats, lons = [], []
    for elem in root.iter():
        if elem.tag == 'latitude':
            lats.append(float(elem.text))
        elif elem.tag == 'longitude':
            lons.append(float(elem.text))
    lat_tie_point = np.mean(lats)
    lon_tie_point = np.mean(lons)
    
    # read LED for SC_vel and SC_height
    sc_vel = 0.0
    x, y, z = 0.0, 0.0, 0.0
    with open(led_file, 'r') as f:
        lines = f.readlines()
        parts = lines[1].split()
        # LED format: year month day x y z vx vy vz
        x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
        vx, vy, vz = float(parts[6]), float(parts[7]), float(parts[8])
        sc_vel = math.sqrt(vx**2 + vy**2 + vz**2)
        sc_h = math.sqrt(x**2 + y**2 + z**2)
        
    a = 6378137.0
    b = 6356752.31
    earth_radius = 6371000.0
    with open(prm_file, 'r') as f:
        lines = f.readlines()
        
    new_lines = []
    for line in lines:
        if line.startswith('SC_vel') or line.startswith('earth_radius') or line.startswith('lon_tie_point') or line.startswith('lat_tie_point') or line.startswith('SC_identity'):
            continue
        new_lines.append(line)
        if 'equatorial_radius' in line:
            a = float(line.split('=')[1].strip())
        elif 'polar_radius' in line:
            b = float(line.split('=')[1].strip())
            
    lat_rad = math.radians(lat_tie_point)
    num = (a**2 * math.cos(lat_rad))**2 + (b**2 * math.sin(lat_rad))**2
    den = (a * math.cos(lat_rad))**2 + (b * math.sin(lat_rad))**2
    earth_radius = math.sqrt(num / den)
    
    sc_height = sc_h
    
    with open(prm_file, 'w') as f:
        for line in new_lines:
            f.write(line)
        f.write(f"lon_tie_point = {lon_tie_point:.6f}\n")
        f.write(f"lat_tie_point = {lat_tie_point:.6f}\n")
        f.write(f"SC_vel = {sc_vel:.6f}\n")
        f.write(f"earth_radius = {earth_radius:.6f}\n")
        f.write(f"SC_height = {sc_height:.6f}\n")
        f.write("SC_identity = 3\n")

if __name__ == "__main__":
    if len(sys.argv) == 4:
        fix_prm(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        # Default fallback for testing
        fix_prm('20240104.PRM', '../raw/bc3-sm-slc-vv-20240104t043852-003394-000033-000d42-01.xml', '20240104.LED')
        fix_prm('20240126.PRM', '../raw/bc3-sm-slc-vv-20240126t043828-003859-000033-000f13-01.xml', '20240126.LED')
