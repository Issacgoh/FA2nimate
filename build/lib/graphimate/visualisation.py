# FA2nimate 

##### Ver:: A0_V1
##### Author(s) : Issac Goh
##### Date : 230813;YYMMDD

### Author notes

### Features to add

def play(filename):
    html = ''
    video = open(filename, 'rb').read()
    src = 'data:video/mp4;base64,' + b64encode(video).decode()
    html += '<video width=1000 controls autoplay loop><source src="%s" type="video/mp4"></video>' % src
    return HTML(html)