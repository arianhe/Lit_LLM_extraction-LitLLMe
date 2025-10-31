import os
import requests

# Zielverzeichnis und Dateiname
model_dir = "./models"
os.makedirs(model_dir, exist_ok=True)
model_path = os.path.join(model_dir, "capybarahermes-2.5-mistral-7b.Q5_K_M.gguf")

# GGUF-Download-Link (direkter Download von HuggingFace)
url = "https://huggingface.co/TheBloke/CapybaraHermes-2.5-Mistral-7B-GGUF/resolve/main/capybarahermes-2.5-mistral-7b.Q5_K_M.gguf"

# Datei herunterladen (nur wenn nicht vorhanden)
if not os.path.exists(model_path):
    print("📥 Lade Modell herunter...")
    response = requests.get(url, stream=True)
    with open(model_path, "wb") as f:
        for chunk in response.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)
    print("✅ Download abgeschlossen:", model_path)
else:
    print("📦 Modell ist bereits vorhanden:", model_path)

# GGUF-Header prüfen
with open(model_path, "rb") as f:
    header = f.read(4)

if header == b'GGUF':
    print("✅ Datei ist gültig: GGUF-Header erkannt.")
else:
    print("❌ Datei ist beschädigt oder kein GGUF-Format!")

