#!/bin/bash
#
# generate_styles.sh - Generate OGC SLD Styling Files for InSAR Products
#
# Part of PowerChina-1 (Spacety) DInSAR System
# Date: 2026-01-04
#
# Output Files:
#   - los_deformation.sld    (Red-White-Blue diverging)
#   - interferogram_phase.sld (Cyclic Rainbow/Jet)
#   - coherence.sld          (Grayscale with transparency)
#   - amplitude.sld          (Linear grayscale)
#

echo "=== Generating OGC SLD Styling Files ==="
echo ""

# Output directory
OUTDIR="${1:-.}"
mkdir -p "$OUTDIR"

#######################################
# 1. LOS Deformation (Red-White-Blue)
#######################################
cat > "$OUTDIR/los_deformation.sld" << 'EOF'
<?xml version="1.0" encoding="UTF-8"?>
<StyledLayerDescriptor version="1.0.0"
    xmlns="http://www.opengis.net/sld"
    xmlns:ogc="http://www.opengis.net/ogc"
    xmlns:xlink="http://www.w3.org/1999/xlink"
    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
    xsi:schemaLocation="http://www.opengis.net/sld http://schemas.opengis.net/sld/1.0.0/StyledLayerDescriptor.xsd">
  <NamedLayer>
    <Name>los_deformation</Name>
    <UserStyle>
      <Title>LOS Deformation (C-band InSAR)</Title>
      <Abstract>Diverging Red-White-Blue color ramp for LOS deformation. Red=Subsidence, Blue=Uplift. Range: -0.05m to +0.05m</Abstract>
      <FeatureTypeStyle>
        <Rule>
          <RasterSymbolizer>
            <ColorMap type="ramp">
              <!-- Subsidence (negative, moving away from satellite) -->
              <ColorMapEntry color="#8B0000" quantity="-0.050" label="-50 mm (Subsidence)" opacity="1.0"/>
              <ColorMapEntry color="#FF0000" quantity="-0.030" label="-30 mm"/>
              <ColorMapEntry color="#FF6666" quantity="-0.015" label="-15 mm"/>
              <ColorMapEntry color="#FFCCCC" quantity="-0.005" label="-5 mm"/>
              <!-- Zero / Stable -->
              <ColorMapEntry color="#FFFFFF" quantity="0.000" label="0 mm (Stable)" opacity="0.8"/>
              <!-- Uplift (positive, moving toward satellite) -->
              <ColorMapEntry color="#CCCCFF" quantity="0.005" label="+5 mm"/>
              <ColorMapEntry color="#6666FF" quantity="0.015" label="+15 mm"/>
              <ColorMapEntry color="#0000FF" quantity="0.030" label="+30 mm"/>
              <ColorMapEntry color="#00008B" quantity="0.050" label="+50 mm (Uplift)" opacity="1.0"/>
            </ColorMap>
          </RasterSymbolizer>
        </Rule>
      </FeatureTypeStyle>
    </UserStyle>
  </NamedLayer>
</StyledLayerDescriptor>
EOF
echo "✅ los_deformation.sld created"

#######################################
# 2. Interferogram Phase (Cyclic Jet)
#######################################
cat > "$OUTDIR/interferogram_phase.sld" << 'EOF'
<?xml version="1.0" encoding="UTF-8"?>
<StyledLayerDescriptor version="1.0.0"
    xmlns="http://www.opengis.net/sld"
    xmlns:ogc="http://www.opengis.net/ogc"
    xmlns:xlink="http://www.w3.org/1999/xlink"
    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
    xsi:schemaLocation="http://www.opengis.net/sld http://schemas.opengis.net/sld/1.0.0/StyledLayerDescriptor.xsd">
  <NamedLayer>
    <Name>interferogram_phase</Name>
    <UserStyle>
      <Title>Interferogram Phase (Wrapped)</Title>
      <Abstract>Cyclic Rainbow (Jet) color ramp for wrapped interferometric phase. Range: -π to +π radians</Abstract>
      <FeatureTypeStyle>
        <Rule>
          <RasterSymbolizer>
            <ColorMap type="ramp">
              <!-- Cyclic: -π to +π (Jet colormap) -->
              <ColorMapEntry color="#00007F" quantity="-3.14159" label="-π" opacity="1.0"/>
              <ColorMapEntry color="#0000FF" quantity="-2.356" label="-3π/4"/>
              <ColorMapEntry color="#00FFFF" quantity="-1.571" label="-π/2"/>
              <ColorMapEntry color="#00FF00" quantity="-0.785" label="-π/4"/>
              <ColorMapEntry color="#FFFF00" quantity="0.000" label="0"/>
              <ColorMapEntry color="#FF7F00" quantity="0.785" label="+π/4"/>
              <ColorMapEntry color="#FF0000" quantity="1.571" label="+π/2"/>
              <ColorMapEntry color="#FF00FF" quantity="2.356" label="+3π/4"/>
              <ColorMapEntry color="#7F0000" quantity="3.14159" label="+π" opacity="1.0"/>
            </ColorMap>
          </RasterSymbolizer>
        </Rule>
      </FeatureTypeStyle>
    </UserStyle>
  </NamedLayer>
</StyledLayerDescriptor>
EOF
echo "✅ interferogram_phase.sld created"

#######################################
# 3. Coherence (Grayscale + Transparency)
#######################################
cat > "$OUTDIR/coherence.sld" << 'EOF'
<?xml version="1.0" encoding="UTF-8"?>
<StyledLayerDescriptor version="1.0.0"
    xmlns="http://www.opengis.net/sld"
    xmlns:ogc="http://www.opengis.net/ogc"
    xmlns:xlink="http://www.w3.org/1999/xlink"
    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
    xsi:schemaLocation="http://www.opengis.net/sld http://schemas.opengis.net/sld/1.0.0/StyledLayerDescriptor.xsd">
  <NamedLayer>
    <Name>coherence</Name>
    <UserStyle>
      <Title>InSAR Coherence</Title>
      <Abstract>Grayscale with transparency for low coherence values. Range: 0.0 (transparent) to 1.0 (white)</Abstract>
      <FeatureTypeStyle>
        <Rule>
          <RasterSymbolizer>
            <ColorMap type="ramp">
              <!-- Low coherence: transparent to indicate unreliable areas -->
              <ColorMapEntry color="#000000" quantity="0.00" label="0.0 (No coherence)" opacity="0.0"/>
              <ColorMapEntry color="#000000" quantity="0.10" label="0.1" opacity="0.3"/>
              <ColorMapEntry color="#1A1A1A" quantity="0.20" label="0.2 (Threshold)" opacity="0.7"/>
              <!-- Medium to high coherence: grayscale -->
              <ColorMapEntry color="#4D4D4D" quantity="0.35" label="0.35" opacity="0.9"/>
              <ColorMapEntry color="#808080" quantity="0.50" label="0.5" opacity="1.0"/>
              <ColorMapEntry color="#B3B3B3" quantity="0.70" label="0.7" opacity="1.0"/>
              <ColorMapEntry color="#E6E6E6" quantity="0.85" label="0.85" opacity="1.0"/>
              <ColorMapEntry color="#FFFFFF" quantity="1.00" label="1.0 (Full coherence)" opacity="1.0"/>
            </ColorMap>
          </RasterSymbolizer>
        </Rule>
      </FeatureTypeStyle>
    </UserStyle>
  </NamedLayer>
</StyledLayerDescriptor>
EOF
echo "✅ coherence.sld created"

#######################################
# 4. Amplitude (Linear Grayscale)
#######################################
cat > "$OUTDIR/amplitude.sld" << 'EOF'
<?xml version="1.0" encoding="UTF-8"?>
<StyledLayerDescriptor version="1.0.0"
    xmlns="http://www.opengis.net/sld"
    xmlns:ogc="http://www.opengis.net/ogc"
    xmlns:xlink="http://www.w3.org/1999/xlink"
    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
    xsi:schemaLocation="http://www.opengis.net/sld http://schemas.opengis.net/sld/1.0.0/StyledLayerDescriptor.xsd">
  <NamedLayer>
    <Name>amplitude</Name>
    <UserStyle>
      <Title>SAR Amplitude</Title>
      <Abstract>Linear grayscale stretch for SAR amplitude imagery. Range: 0 (black) to 255 (white)</Abstract>
      <FeatureTypeStyle>
        <Rule>
          <RasterSymbolizer>
            <ColorMap type="ramp">
              <!-- Linear grayscale: 0-255 -->
              <ColorMapEntry color="#000000" quantity="0" label="0 (Min)" opacity="1.0"/>
              <ColorMapEntry color="#1A1A1A" quantity="25" label="25"/>
              <ColorMapEntry color="#333333" quantity="50" label="50"/>
              <ColorMapEntry color="#4D4D4D" quantity="75" label="75"/>
              <ColorMapEntry color="#666666" quantity="100" label="100"/>
              <ColorMapEntry color="#808080" quantity="125" label="125"/>
              <ColorMapEntry color="#999999" quantity="150" label="150"/>
              <ColorMapEntry color="#B3B3B3" quantity="175" label="175"/>
              <ColorMapEntry color="#CCCCCC" quantity="200" label="200"/>
              <ColorMapEntry color="#E6E6E6" quantity="225" label="225"/>
              <ColorMapEntry color="#FFFFFF" quantity="255" label="255 (Max)" opacity="1.0"/>
            </ColorMap>
          </RasterSymbolizer>
        </Rule>
      </FeatureTypeStyle>
    </UserStyle>
  </NamedLayer>
</StyledLayerDescriptor>
EOF
echo "✅ amplitude.sld created"

#######################################
# Summary
#######################################
echo ""
echo "============================================="
echo "  SLD Styling Files Generated Successfully"
echo "============================================="
echo ""
echo "Output Directory: $OUTDIR"
echo ""
echo "Files Created:"
ls -la "$OUTDIR"/*.sld 2>/dev/null
echo ""
echo "Usage with GeoServer/MapServer:"
echo "  1. Upload .sld files to your map server"
echo "  2. Associate each style with corresponding layer:"
echo "     - los_deformation.sld   → LOS deformation GeoTIFF"
echo "     - interferogram_phase.sld → Phase/phasefilt GeoTIFF"
echo "     - coherence.sld         → Coherence GeoTIFF"
echo "     - amplitude.sld         → Amplitude GeoTIFF"
echo ""
