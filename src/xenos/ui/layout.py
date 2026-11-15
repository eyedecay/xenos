"""Shared Xenos page styling and lightweight display helpers."""

def render_page_chrome(st):
    st.set_page_config(page_title="Xenos | Planetary analog explorer", page_icon="◉", layout="wide")
    st.markdown("""<style>
@import url('https://fonts.googleapis.com/css2?family=DM+Mono:wght@400;500&family=DM+Sans:wght@400;500;600;700&display=swap');
html, body, [class*="css"] {font-family:'DM Sans',sans-serif}
.stApp {background:#0b1013;color:#e7e9e7}
[data-testid="stSidebar"] {background:#11191d;border-right:1px solid #263238}
.eyebrow {color:#a5b9b7;font:500 11px 'DM Mono',monospace;letter-spacing:1.6px;text-transform:uppercase}
.hero {padding:22px 0 16px;border-bottom:1px solid #263238;margin-bottom:18px}
.hero h1 {font-size:34px;letter-spacing:-1.3px;margin:5px 0 4px;color:#f4f2ea}
.hero p {color:#9eaaa9;margin:0;font-size:14px}
.metric {background:#131c20;border:1px solid #28363b;border-radius:10px;padding:15px 17px;min-height:91px}
.metric .label {font:11px 'DM Mono',monospace;color:#95a6a5;text-transform:uppercase;letter-spacing:.8px}
.metric .value {font-size:25px;font-weight:600;color:#f4f2ea;margin-top:5px}
.metric .unit {font:12px 'DM Mono',monospace;color:#91a19f}
.panel {background:#11191d;border:1px solid #263238;border-radius:12px;padding:17px 18px;margin-bottom:12px}
.panel-title {font-size:16px;font-weight:600;color:#f0efe8;margin-bottom:3px}
.panel-sub {font-size:12px;color:#8f9e9e;margin-bottom:14px}
.tag {display:inline-block;border:1px solid #324247;border-radius:20px;padding:4px 9px;color:#abc0bd;font:10px 'DM Mono',monospace;margin-right:5px}
div.stButton>button {background:#c5f36b;color:#101811;border:0;border-radius:7px;font-weight:700;padding:10px 16px}
div.stButton>button:hover {background:#d7ff8a;color:#101811}
.smallnote {font-size:11px;line-height:1.55;color:#889896}
</style>""", unsafe_allow_html=True)
    st.markdown('<div class="hero"><div class="eyebrow">NASA Space Apps · planetary field intelligence</div><h1>Xenos <span style="color:#c5f36b">/</span></h1><p>Discover and compare Earth sites for planetary terrain and geology field tests.</p></div>', unsafe_allow_html=True)

def metric(st, label, value, unit):
    st.markdown(f'<div class="metric"><div class="label">{label}</div><div class="value">{value} <span class="unit">{unit}</span></div></div>', unsafe_allow_html=True)

def format_coordinates(lat, lon, precision=3):
    north_south = "N" if lat >= 0 else "S"
    east_west = "E" if lon >= 0 else "W"
    return f"{abs(lat):.{precision}f}°{north_south}, {abs(lon):.{precision}f}°{east_west}"
