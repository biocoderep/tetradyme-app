import json
from fastapi import FastAPI, UploadFile, Form, File, HTTPException
from fastapi.responses import JSONResponse
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from core import process_primers

app = FastAPI(title="Tetradyme API")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Mount static files to serve the frontend
app.mount("/static", StaticFiles(directory="static"), name="static")

@app.post("/api/design")
async def design_primers(
    fasta_file: UploadFile = File(...),
    csv_file: UploadFile = File(...),
    na: float = Form(50.0),
    mg: float = Form(0.0),
    dntps: float = Form(0.0),
    tm_tolerance: float = Form(1.0),
    of_primer: str = Form(""),
    or_primer: str = Form("")
):
    if not fasta_file.filename.endswith((".fasta", ".fa")):
        raise HTTPException(status_code=400, detail="FASTA file must have .fasta or .fa extension")
    if not csv_file.filename.endswith(".csv"):
        raise HTTPException(status_code=400, detail="CSV file must have .csv extension")

    fasta_content = (await fasta_file.read()).decode("utf-8")
    csv_content = (await csv_file.read()).decode("utf-8")

    params = {
        "na": na,
        "mg": mg,
        "dntps": dntps,
        "tm_tolerance": tm_tolerance,
        "of_primer": of_primer,
        "or_primer": or_primer
    }

    try:
        result = process_primers(fasta_content, csv_content, params)
        return JSONResponse(content=result)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@app.get("/")
async def root():
    # Redirect root to index.html in static
    from fastapi.responses import RedirectResponse
    return RedirectResponse(url="/static/index.html")
