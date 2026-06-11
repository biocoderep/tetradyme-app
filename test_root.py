import urllib.request
try:
    with urllib.request.urlopen('http://127.0.0.1:8000/') as response:
        print("Root Response Code:", response.getcode())
        print("Root Response URL (after redirect):", response.geturl())
except Exception as e:
    print("Error:", e)
