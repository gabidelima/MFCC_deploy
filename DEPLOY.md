# Deploy no Streamlit Community Cloud

## Passo a passo (5 minutos, gratuito)

### 1. Suba o repositório no GitHub

Estrutura obrigatória:
```
mfcc-automation/
├── app.py                    ← ponto de entrada do Streamlit
├── mfcc_engine.py            ← motor Python puro
├── generate_database.py
├── fragment_by_distance.py
├── requirements.txt
└── packages.txt
```

> Os scripts `MFCC_CPCM_*.sh` **não são mais necessários** para o deploy em nuvem.
> O `mfcc_engine.py` replaca toda a lógica bash em Python puro.

### 2. Crie uma conta em share.streamlit.io

Acesse: https://share.streamlit.io  
Faça login com sua conta GitHub.

### 3. Clique em "New app"

Preencha:
- **Repository**: `seu-usuario/mfcc-automation`
- **Branch**: `main`
- **Main file path**: `app.py`

Clique em **Deploy**.

### 4. Aguarde ~2 minutos

O Streamlit Cloud instala as dependências automaticamente a partir do `requirements.txt`.

### 5. URL pública gerada

```
https://seu-usuario-mfcc-automation-app-main.streamlit.app
```

Compartilhe essa URL com qualquer pessoa — sem login, sem instalação.

---

## Limites do plano gratuito

| Recurso         | Limite                        |
|----------------|-------------------------------|
| Apps públicos   | Ilimitados                    |
| RAM             | 1 GB por app                  |
| CPU             | Compartilhada                 |
| Armazenamento   | Temporário (sem persistência) |
| Sleep automático| App dorme após inatividade    |

> Para sistemas grandes (>300 resíduos, >12 Å), o processamento pode
> ser lento na nuvem gratuita. Para uso intensivo, considere um VPS.

---

## Rodar localmente (alternativa)

```bash
git clone https://github.com/seu-usuario/mfcc-automation.git
cd mfcc-automation
pip install -r requirements.txt
streamlit run app.py
```

Abre em: http://localhost:8501
