mkdir -p ~/.streamlit/
echo "[general]  
email = \"me@gmail.com\""  > ~/.streamlit/credentials.toml
echo "[server]
headless = true
port = $PORT
enableCORS = true"  >> ~/.streamlit/config.toml
