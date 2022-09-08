import io
import streamlit
import numpy
import base64

# import scipy
# import scipy.optimize

def render_svg(svg, col):
    """Renders the given svg string."""
    b64 = base64.b64encode(svg.encode('utf-8')).decode("utf-8")
    # Add some CSS on top
    center = True
    css_justify = "center" if center else "left"
    css = '<p style="text-align:center; display: flex; justify-content: {};">'.format(css_justify)
    html = r'{}<img src="data:image/svg+xml;base64,{}"/>'.format(
        css, b64)
        
    # html = r'%s' % b64
    col.write(html, unsafe_allow_html=True)
    return



def svg_circle(x,y,rx):
    ls_out = []
    return "\n".join(ls_out)


def from_coord_to_frame_1d(coord_limit, frame_limit, coord):
    coord_min = coord_limit[0]
    delta_coord = coord_limit[1] - coord_min

    frame_min = frame_limit[0]
    delta_frame = frame_limit[1] - frame_min

    frame = (coord - coord_min) * delta_frame / delta_coord
    return frame

def from_coord_to_frame_2d(coord_limit, frame_limit, coord):
    coord_limit_x = (coord_limit[0], coord_limit[2])
    coord_limit_y = (coord_limit[1], coord_limit[3])
    frame_limit_x = (frame_limit[0], frame_limit[2])
    frame_limit_y = (frame_limit[1], frame_limit[3])
    coord_x = coord[0]
    coord_y = -coord[1]
    frame_x = from_coord_to_frame_1d(coord_limit_x, frame_limit_x, coord_x)
    frame_y = from_coord_to_frame_1d(coord_limit_y, frame_limit_y, coord_y)
    frame = numpy.stack([frame_x, frame_y], axis=0)
    return frame


def svg_head(frame_limit):
    ls_out = []
    ls_out.append(f"<svg width=\"10cm\" height=\"10cm\" viewBox=\"{frame_limit[0]:} {frame_limit[1]:} {frame_limit[2]:} {frame_limit[3]:}\"")
    s_tail = """xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" xml:space="preserve">
<title>Example of path</title>
<defs>
<marker id='head' orient="auto"
        markerWidth='3' markerHeight='4'
        refX='0.1' refY='2'>
<path d='M0,0 V4 L2,2 Z' fill="black"/>
</marker>

</defs>
"""
    ls_out.append(s_tail)
    return "\n".join(ls_out)

def svg_tail():
    s_out = "</svg>"
    return s_out

def svg_circle(frame_xy, **d_arg):
    d_arg_keys = d_arg.keys()
    color = "#000000"
    if "color" in d_arg_keys:
        color = d_arg["color"]


    ls_circle = []
    for x, y in zip(frame_xy[0], frame_xy[1]):
        s_circle = f"<circle id=\"point_b\" r=\"2\" cx=\"{x:.2f}\" cy=\"{y:.2f}\" opacity=\"0.2\" fill=\"none\" stroke=\"{color:}\">\n</circle>\n"
        ls_circle.append(s_circle)
    return "\n".join(ls_circle)


def svg_animated_line_from_center(frame_xy, frame_center, **d_arg):
    d_arg_keys = d_arg.keys()
    color = "#000000"
    if "color" in d_arg_keys:
        color = d_arg["color"]

    ls_x, ls_y = [], []
    for x, y in zip(frame_xy[0], frame_xy[1]):
        ls_x.append(f"{x:.2f}")
        ls_y.append(f"{y:.2f}")
    s_x = ";".join(ls_x)
    s_y = ";".join(ls_y)

    ls_out = []
    ls_out.append(f"<line x1=\"{frame_center[0]:.2f}\" y1=\"{frame_center[1]:.2f}\" x2=\"{frame_xy[0,0]:.2f}\" y2=\"{frame_xy[1,0]:.2f}\" stroke=\"{color:}\" marker-end='url(#head)'>")
    ls_out.append(f"<animate attributeName=\"x2\" values=\"{s_x:}\" dur=\"9s\" repeatCount=\"indefinite\" />")
    ls_out.append(f"<animate attributeName=\"y2\" values=\"{s_y:}\" dur=\"9s\" repeatCount=\"indefinite\" />")
    ls_out.append("</line>")
    return "\n".join(ls_out)

def svg_animated_text(frame_xy, s_text: str, **d_arg):
    d_arg_keys = d_arg.keys()
    color = "#000000"
    if "color" in d_arg_keys:
        color = d_arg["color"]

    ls_x, ls_y = [], []
    for x, y in zip(frame_xy[0], frame_xy[1]):
        ls_x.append(f"{x:.2f}")
        ls_y.append(f"{y:.2f}")
    s_x = ";".join(ls_x)
    s_y = ";".join(ls_y)

    ls_out = []
    ls_out.append(f"<text x=\"{frame_xy[0,0]:.2f}\" y=\"{frame_xy[1,0]:.2f}\" font-size=\"12\" fill=\"{color:}\">{s_text}")
    ls_out.append(f"<animate attributeName=\"x\" values=\"{s_x:}\" dur=\"9s\" repeatCount=\"indefinite\" />")
    ls_out.append(f"<animate attributeName=\"y\" values=\"{s_y:}\" dur=\"9s\" repeatCount=\"indefinite\" />")
    ls_out.append("</text>")

    return "\n".join(ls_out)

def animated_block(frame_limit, coord_limit, coord_b, label_b, **d_arg):
    # simplified expression for the center
    frame_center = (0.5 * frame_limit[0]+ 0.5 * frame_limit[2], 0.5 * frame_limit[1]+ 0.5 * frame_limit[3])
    frame_b = from_coord_to_frame_2d(coord_limit, frame_limit, coord_b)
    svg_circle_b =  svg_circle(frame_b, **d_arg)
    svg_line = svg_animated_line_from_center(frame_b, frame_center, **d_arg)
    svg_text = svg_animated_text(frame_b, label_b, **d_arg)
    ls_out = [svg_circle_b, svg_line, svg_text]
    return "\n".join(ls_out)

def coord_limit(*argv):
    border = 1.1
    l_xy_max = []
    for coord_xy in argv:
        abs_coord_xy = numpy.abs(coord_xy)
        l_xy_max.append(numpy.max(abs_coord_xy, axis=1))
    np_xy_max = numpy.stack(l_xy_max, axis=1)
    xy_max = numpy.max(np_xy_max)*border
    coord_limit = (-xy_max, -xy_max, xy_max, xy_max)
    return coord_limit





def calc_xyz_by_theta_phi(theta, phi):
    h_x = numpy.cos(phi) * numpy.sin(theta)
    h_y = numpy.sin(phi) * numpy.sin(theta)
    h_z = numpy.cos(theta) * numpy.ones_like(phi)
    return h_x, h_y, h_z


def calc_energy_uniaxial_symmetry(theta, phi, k_0, k_1):
    h_x, h_y, h_z = calc_xyz_by_theta_phi(theta, phi)
    energy = 0.5 * (k_1 * (numpy.square(h_x) + numpy.square(h_y)) + k_0 * numpy.square(h_z))
    return energy


def calc_moment_uniaxial_symmetry(theta, phi, k_0, k_1):
    h_x, h_y, h_z = calc_xyz_by_theta_phi(theta, phi)
    m_x = k_1 * h_x
    m_y = k_1 * h_y
    m_z = k_0 * h_z
    m = numpy.stack([m_x, m_y, m_z], axis=0)
    return m


def calc_energy_rhombohedral_symmetry(theta, phi, k_0, k_1, k_2, k_3, k_4):
    h_x, h_y, h_z = calc_xyz_by_theta_phi(theta, phi)
    h_x_sq, h_y_sq, h_z_sq = numpy.square(h_x), numpy.square(h_y), numpy.square(h_z)

    energy = 0.5 * (
        k_1 * (numpy.square(h_x) + numpy.square(h_y)) + 
        k_0 * numpy.square(h_z)) + (
            k_2 * (h_x_sq - 3*h_y_sq) * h_x +
            k_3 * (h_y_sq - 3*h_x_sq) * h_y + 
            k_4 * (5 * h_z_sq - 3) * h_z
        )
    return energy


def calc_moment_rhombohedral_symmetry(theta, phi, k_0, k_1, k_2, k_3, k_4):
    h_x, h_y, h_z = calc_xyz_by_theta_phi(theta, phi)
    h_x_sq, h_y_sq, h_z_sq = numpy.square(h_x), numpy.square(h_y), numpy.square(h_z)

    m_x = 2 * 0.5 * k_1 * h_x + 3 * k_2 * (h_x_sq - h_y_sq) - 6 * k_3 * h_x * h_y
    m_y = 2 * 0.5 * k_1 * h_y - 6 * k_2 * h_x * h_y + 3 * k_3 * (h_y_sq - h_x_sq)
    m_z = 2 * 0.5 * k_0 * h_z + k_4 * (5 * h_z_sq + 10 * h_z - 3) 
    m = numpy.stack([m_x, m_y, m_z], axis=0)
    return m

def calc_moment_tetragonal_symmetry(theta, phi, k_0, k_1, k_2, k_3, k_4):
    h_x, h_y, h_z = calc_xyz_by_theta_phi(theta, phi)
    h_x_sq, h_y_sq, h_z_sq = numpy.square(h_x), numpy.square(h_y), numpy.square(h_z)

    m_x = 2 * 0.5 * k_1 * h_x + 4 * k_2 * h_x_sq * h_x + 2 * k_3 * h_x * h_y_sq 
    m_y = 2 * 0.5 * k_1 * h_y + 4 * k_2 * h_y_sq * h_y + 2 * k_3 * h_x_sq * h_y
    m_z = 2 * 0.5 * k_0 * h_z + 4*k_4 * h_z_sq * h_z 
    m = numpy.stack([m_x, m_y, m_z], axis=0)
    return m

def calc_moment_ortogonal_symmetry(theta, phi, k_0, k_1, k_2):
    h_x, h_y, h_z = calc_xyz_by_theta_phi(theta, phi)
    h_x_sq, h_y_sq, h_z_sq = numpy.square(h_x), numpy.square(h_y), numpy.square(h_z)

    m_x = 2 * 0.5 * k_1 * h_x  
    m_y = 2 * 0.5 * k_2 * h_y 
    m_z = 2 * 0.5 * k_0 * h_z  
    m = numpy.stack([m_x, m_y, m_z], axis=0)
    return m


def calc_fields_and_moments(dict_k: dict, n_points: int = 90):
    theta, phi = 0.5*numpy.pi, 0.
    np_angle_phi = numpy.linspace(0., 2.*numpy.pi, n_points, endpoint=False)
    np_angle_theta =  - np_angle_phi
    b_theta = numpy.stack(calc_xyz_by_theta_phi(np_angle_theta, phi), axis=0)
    b_phi = numpy.stack(calc_xyz_by_theta_phi(theta, np_angle_phi), axis=0)
    
    
    model = dict_k["model"]
    if (model == "uniaxial"):
        k_0 = dict_k["k_0"]
        k_1 = dict_k["k_1"]
        m_over_theta = calc_moment_uniaxial_symmetry(np_angle_theta, phi, k_0, k_1)
        energy_over_theta = calc_moment_uniaxial_symmetry(np_angle_theta, phi, k_0, k_1)

        m_over_phi = calc_moment_uniaxial_symmetry(theta, np_angle_phi, k_0, k_1)
        energy_over_phi = calc_moment_uniaxial_symmetry(theta, np_angle_phi, k_0, k_1)
    elif (model == "rhombohedral"):
        k_0 = dict_k["k_0"]
        k_1 = dict_k["k_1"]
        k_2 = dict_k["k_2"]
        k_3 = dict_k["k_3"]
        k_4 = dict_k["k_4"]

        m_over_theta = calc_moment_rhombohedral_symmetry(np_angle_theta, phi, k_0, k_1, k_2, k_3, k_4)
        energy_over_theta = calc_moment_rhombohedral_symmetry(np_angle_theta, phi, k_0, k_1, k_2, k_3, k_4)

        m_over_phi = calc_moment_rhombohedral_symmetry(theta, np_angle_phi, k_0, k_1, k_2, k_3, k_4)
        energy_over_phi = calc_moment_rhombohedral_symmetry(theta, np_angle_phi, k_0, k_1, k_2, k_3, k_4)
    elif (model == "tetragonal"):
        k_0 = dict_k["k_0"]
        k_1 = dict_k["k_1"]
        k_2 = dict_k["k_2"]
        k_3 = dict_k["k_3"]
        k_4 = dict_k["k_4"]

        m_over_theta = calc_moment_tetragonal_symmetry(np_angle_theta, phi, k_0, k_1, k_2, k_3, k_4)
        energy_over_theta = calc_moment_tetragonal_symmetry(np_angle_theta, phi, k_0, k_1, k_2, k_3, k_4)

        m_over_phi = calc_moment_tetragonal_symmetry(theta, np_angle_phi, k_0, k_1, k_2, k_3, k_4)
        energy_over_phi = calc_moment_tetragonal_symmetry(theta, np_angle_phi, k_0, k_1, k_2, k_3, k_4)
    elif (model == "ortogonal"):
        k_0 = dict_k["k_0"]
        k_1 = dict_k["k_1"]
        k_2 = dict_k["k_2"]

        m_over_theta = calc_moment_ortogonal_symmetry(np_angle_theta, phi, k_0, k_1, k_2)
        energy_over_theta = calc_moment_ortogonal_symmetry(np_angle_theta, phi, k_0, k_1, k_2)

        m_over_phi = calc_moment_ortogonal_symmetry(theta, np_angle_phi, k_0, k_1, k_2)
        energy_over_phi = calc_moment_ortogonal_symmetry(theta, np_angle_phi, k_0, k_1, k_2)
    else:
        raise KeyError("defined models: 'uniaxial' or 'rhombohedral'")

    return b_theta, m_over_theta, energy_over_theta, b_phi, m_over_phi, energy_over_phi


def get_energy_uniaxial():
    str_out = r"E = \frac{K_0  B_z^2 + K_1 (B_x^2 + B_y^2)}{2}"
    return str_out


def get_energy_rhombohedral():
    str_out = r"E = \frac{K_0  B_z^2 + K_1 (B_x^2 + B_y^2)}{2}  + K_2 ( B_x^2 - B_y^2 ) B_x + K_3 ( B_y^2 - B_x^2 ) B_y + K_4 ( 5 B_z^2 - 3 ) B_z"
    return str_out


def get_energy_tetragonal():
    str_out = r"E = \frac{K_0  B_z^2 + K_1 (B_x^2 + B_y^2)}{2}  + K_2 ( B_x^4 + B_y^4 ) + K_3 B_x^2 B_y^2  + K_4 B_z^4"
    return str_out

def get_energy_ortogonal():
    str_out = r"E = \frac{K_0  B_z^2 + K_1 B_x^2 + K_2 B_y^2}{2}"
    return str_out


streamlit.markdown("Calculate magnetic moment $\mathbf{M}$ of one ion in the following model: ")

d_params = {}

d_params["model"] = streamlit.selectbox("", ["uniaxial", "rhombohedral", "tetragonal", "ortogonal"], index=0)

if d_params["model"] == "uniaxial":
    streamlit.latex(get_energy_uniaxial())
elif d_params["model"] == "rhombohedral":
    streamlit.latex(get_energy_rhombohedral())
elif d_params["model"] == "tetragonal":
    streamlit.latex(get_energy_tetragonal())
elif d_params["model"] == "ortogonal":
    streamlit.latex(get_energy_ortogonal())

streamlit.markdown("Parameters of the model: ")

col_k0, col_k1, col_k2 = streamlit.columns(3)
d_params["k_0"] = col_k0.number_input(f"K0:", -10., 10., value=1., step=0.1)
d_params["k_1"] = col_k1.number_input(f"K1:", -10., 10., value=0., step=0.1)
d_params["k_2"] = col_k2.number_input(f"K2:", -10., 10., value=0., step=0.1)
d_params["k_3"] = col_k0.number_input(f"K3:", -10., 10., value=0., step=0.1)
d_params["k_4"] = col_k1.number_input(f"K4:", -10., 10., value=0., step=0.1)
d_params["k_5"] = col_k2.number_input(f"K5:", -10., 10., value=0., step=0.1)


# theta = streamlit.number_input("Theta", 0., 360., value=0., step=15.) * numpy.pi/180.
# phi = streamlit.number_input("Phi", 0., 360., value=0., step=15.) * numpy.pi/180.

b_moments = streamlit.button("Calculate")
if b_moments:
    with streamlit.spinner("Please wait..."):
        n_points = 180 
        b_theta, m_over_theta, energy_over_theta, b_phi, m_over_phi, energy_over_phi = \
            calc_fields_and_moments(d_params, n_points = n_points)

    col_1, col_2 = streamlit.columns(2)

    np_b_x, np_b_z = b_phi[0], b_phi[1] 
    np_m_x, np_m_z = m_over_phi[0], m_over_phi[1]

    coord_b = numpy.stack([np_b_x, np_b_z], axis=0)
    coord_m = numpy.stack([np_m_x, np_m_z], axis=0)

    coord_l = coord_limit(coord_b, coord_m)
    frame_limit = (0, 0, 150, 150)
    label_b = "B_xy"
    label_m = "m_xy"
    ls_svg = []
    ls_svg.append(svg_head(frame_limit))

    svg_block_b = animated_block(frame_limit, coord_l, coord_b, label_b, color="#ff0000")
    ls_svg.append(svg_block_b)

    svg_block_m = animated_block(frame_limit, coord_l, coord_m, label_m, color="#0000ff")
    ls_svg.append(svg_block_m)

    ls_svg.append(svg_tail())

    render_svg("\n".join(ls_svg), col_1)



    np_b_x, np_b_z = b_theta[0], b_theta[2] 
    np_m_x, np_m_z = m_over_theta[0], m_over_theta[2]

    coord_b = numpy.stack([np_b_x, np_b_z], axis=0)
    coord_m = numpy.stack([np_m_x, np_m_z], axis=0)

    coord_l = coord_limit(coord_b, coord_m)
    frame_limit = (0, 0, 150, 150)
    label_b = "B_xz"
    label_m = "m_xz"
    ls_svg = []
    ls_svg.append(svg_head(frame_limit))

    svg_block_b = animated_block(frame_limit, coord_l, coord_b, label_b, color="#ff0000")
    ls_svg.append(svg_block_b)

    svg_block_m = animated_block(frame_limit, coord_l, coord_m, label_m, color="#0000ff")
    ls_svg.append(svg_block_m)

    ls_svg.append(svg_tail())

    render_svg("\n".join(ls_svg), col_2)

streamlit.markdown("""
Please visit [the site][link_streamlit] to estimate the effect of magnetic field on interacting ions.

[link_streamlit]: https://ikibalin-moments-main-0doyyd.streamlitapp.com/
""")