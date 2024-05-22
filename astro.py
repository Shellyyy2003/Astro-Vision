from astropy.table import Table
import requests
from PIL import Image
from io import BytesIO
from flask import Flask, render_template, request, send_file, jsonify
from io import BytesIO
import numpy
from astropy.table import Table
import requests
from PIL import Image
from flask import send_file
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from io import BytesIO
from matplotlib.colors import Normalize
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time

latest_csv_data = None

def load_default_csv():
    global latest_csv_data
    default_csv_path = 'candidatesNotInBaseLG_visualMark.csv'
    latest_csv_data = pd.read_csv(default_csv_path).head(100)

def getimages(ra,dec,filters="grizy"):
    
    """Query ps1filenames.py service to get a list of images
    
    ra, dec = position in degrees
    size = image size in pixels (0.25 arcsec/pixel)
    filters = string with filters to include
    Returns a table with the results
    """
    
    service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = f"{service}?ra={ra}&dec={dec}&filters={filters}"
    table = Table.read(url, format='ascii')
    return table


def geturl(ra, dec, size=240, output_size=None, filters="grizy", format="jpg", color=False):
    
    """Get URL for images in the table
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png" or "fits")
    color = if True, creates a color image (only for jpg or png format).
            Default is return a list of URLs for single-filter grayscale images.
    Returns a string with the URL
    """
    
    if color and format == "fits":
        raise ValueError("color images are available only for jpg or png formats")
    if format not in ("jpg","png","fits"):
        raise ValueError("format must be one of jpg, png, fits")
    table = getimages(ra,dec,filters=filters)
    url = (f"https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
           f"ra={ra}&dec={dec}&size={size}&format={format}")
    if output_size:
        url = url + "&output_size={}".format(output_size)
    # sort filters from red to blue
    flist = ["grizy".find(x) for x in table['filter']]
    table = table[numpy.argsort(flist)]
    if color:
        if len(table) > 3:
            # pick 3 filters
            table = table[[0,len(table)//2,len(table)-1]]
        for i, param in enumerate(["red","green","blue"]):
            url = url + "&{}={}".format(param,table['filename'][i])
    else:
        urlbase = url + "&red="
        url = []
        for filename in table['filename']:
            url.append(urlbase+filename)
    return url


def getcolorim(ra, dec, size=240, output_size=None, filters="grizy", format="jpg"):
    
    """Get color image at a sky position
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png")
    Returns the image
    """
    
    if format not in ("jpg","png"):
        raise ValueError("format must be jpg or png")
    url = geturl(ra,dec,size=size,filters=filters,output_size=output_size,format=format,color=True)
    r = requests.get(url)
    im = Image.open(BytesIO(r.content))
    return im


def getgrayim(ra, dec, size=240, output_size=None, filter="g", format="jpg"):
    
    """Get grayscale image at a sky position
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filter = string with filter to extract (one of grizy)
    format = data format (options are "jpg", "png")
    Returns the image
    """
    
    if format not in ("jpg","png"):
        raise ValueError("format must be jpg or png")
    if filter not in list("grizy"):
        raise ValueError("filter must be one of grizy")
    url = geturl(ra,dec,size=size,filters=filter,output_size=output_size,format=format)
    r = requests.get(url[0])
    im = Image.open(BytesIO(r.content))
    return im


load_default_csv()
app = Flask(__name__)

FILTER_COLORS = {
    "g": "#78A5A3",
    "r": "#CE5A57",
    "i": "#FF496C",
    "z": "#785A46",
    "y": "#F9E076"
}

# Route for the main page which displays the HTML form
@app.route('/')
def display_form():
    return render_template('index.html')

# Route to handle the image URL fetching based on user input
@app.route('/geturl', methods=['POST'])
def get_image_url():
    data = request.get_json()
    
    print("In get_image_url, data is", data)
    if isinstance(data['ra'], dict):
        ra = float(data['ra'].get('RA') or data['ra'].get('ra'))
    else:
        ra = float(data['ra'])
    
    if isinstance(data['dec'], dict):
        dec = float(data['dec'].get('Dec') or data['dec'].get('DEC') or data['dec'].get('dec'))
    else:
        dec = float(data['dec'])
    size = int(data['size']) if 'size' in data and data['size'] else 240
    
    url_list = geturl(ra, dec, size, filters="grizy", color=False)
    images = [{"url": url, "filter": filter, "color": FILTER_COLORS[filter]} for url, filter in zip(url_list, "grizy")]
    return jsonify(images=images)


@app.route('/get_star_chart', methods=['POST'])
def get_star_chart():
    try:
        data = request.get_json()
        print("data.get is: ",data.get('ra', 0))
        if isinstance(data['ra'], dict):
            ra = float(data['ra'].get('RA') or data['ra'].get('ra'))
        else:
            ra = float(data['ra'])
        
        if isinstance(data['dec'], dict):
            dec = float(data['dec'].get('Dec') or data['dec'].get('DEC') or data['dec'].get('dec'))
        else:
            dec = float(data['dec'])
        highlight_index = data.get('highlightIndex', -1)
        start_time = time.time()

        global latest_csv_data

        # 创建图表并绘制 CSV 中的点
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(111, projection="aitoff")

        for index, row in latest_csv_data.iterrows():
            if time.time() - start_time > 5:  # Check if 10 seconds have passed
                break
            c = SkyCoord(ra=row['RA'] * u.deg, dec=row['Dec'] * u.deg, frame='icrs').galactic
            color = 'red' if index == highlight_index else 'gray'
            ax.scatter(c.l.wrap_at(180 * u.deg).radian, c.b.radian, s=10, color=color)

        # 如果用户输入了数据，绘制这个点
        if ra != 0 or dec != 0:
            c_user = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs').galactic
            ax.scatter(c_user.l.wrap_at(180 * u.deg).radian, c_user.b.radian, s=50, color='blue')

        ax.set_xlabel('Galactic Longitude [degrees]')
        ax.set_ylabel('Galactic Latitude [degrees]')
        ax.grid(True)

        buf = BytesIO()
        plt.savefig(buf, format='png')
        buf.seek(0)
        return send_file(buf, mimetype='image/png')
    except TypeError as e:
        print(f"An error occurred: {e}")
        return jsonify(error=str(e)), 400


@app.route('/get_unchanged_star_chart', methods=['POST'])
def get_unchanged_star_chart():
    data = request.get_json()
    if isinstance(data['ra'], dict):
        ra = float(data['ra'].get('RA') or data['ra'].get('ra'))
    else:
        ra = float(data['ra'])
    
    if isinstance(data['dec'], dict):
        dec = float(data['dec'].get('Dec') or data['dec'].get('DEC') or data['dec'].get('dec'))
    else:
        dec = float(data['dec'])
    highlight_index = data.get('highlightIndex', -1)
    start_time = time.time()

    global latest_csv_data

    # 创建图表并绘制 CSV 中的点
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111, projection="aitoff")

    for index, row in latest_csv_data.iterrows():
        if time.time() - start_time > 5:  # Check if 10 seconds have passed
                break
        c = SkyCoord(ra=row['RA'] * u.deg, dec=row['Dec'] * u.deg, frame='icrs').galactic
        color = 'red' if index == highlight_index else 'gray'
        ax.scatter(c.l.wrap_at(180 * u.deg).radian, c.b.radian, s=10, color=color)

    # 如果用户输入了数据，绘制这个点
    if ra != 0 or dec != 0:
        c_user = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs')
        ax.scatter(c_user.ra.wrap_at(180 * u.deg).radian, c_user.dec.radian, s=50, color='blue')

    # 设置标签和网格
    ax.set_xlabel('Right Ascension [degrees]')
    ax.set_ylabel('Declination [degrees]')
    ax.grid(True)

    buf = BytesIO()
    plt.savefig(buf, format='png')
    buf.seek(0)
    return send_file(buf, mimetype='image/png')


@app.route('/get_csv_data', methods=['GET'])
def get_csv_data():
    global latest_csv_data
    columns = latest_csv_data.columns.tolist()
    return jsonify(data=latest_csv_data.to_json(orient='records'), columns=columns)

def create_initial_star_chart():
    global latest_csv_data

    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111, projection="aitoff")

    # 绘制 CSV 文件中的点
    for index, row in latest_csv_data.iterrows():
        c = SkyCoord(ra=row['RA'] * u.deg, dec=row['Dec'] * u.deg, frame='icrs')
        galactic = c.galactic
        wrapped_l = galactic.l.wrap_at(180 * u.deg).radian
        b = galactic.b.radian
        ax.scatter(wrapped_l, b, s=10, color='red')

    # 设置标签和网格
    ax.set_xlabel('Galactic Longitude [degrees]')
    ax.set_ylabel('Galactic Latitude [degrees]')
    ax.grid(True)

    # 保存图像
    buf = BytesIO()
    plt.savefig(buf, format='png')
    buf.seek(0)
    return buf.getvalue()

@app.route('/upload_csv', methods=['POST'])
def upload_csv():
    if 'file' not in request.files:
        print("File not in request.files")
        return jsonify(success=False, message="No file part")
    
    file = request.files['file']
    if file.filename == '':
        print("File name is empty")
        return jsonify(success=False, message="No selected file")
    
    if file:
        try:
            global latest_csv_data
            latest_csv_data = pd.read_csv(file)
            return jsonify(success=True, message="File uploaded successfully")
        except Exception as e:
            return jsonify(success=False, message=str(e))
    
    return jsonify(success=False, message="Invalid file format")


if __name__ == "__main__":
    app.run(debug=True, threaded=False)
