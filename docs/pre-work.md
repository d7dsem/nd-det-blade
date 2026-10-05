# bladeRF + Python на Windows 10 — шпаргалка

Реліз bladeRF: **2025.10** (libbladeRF 2.6.0, FX3 2.6.0, FPGA 0.16.0, Python bindings 1.4.0). Перевірено на: Win10 22H2, Python 3.13.9 x64, bladeRF 2.0 micro.

---

## 1. Перевірка оточення (лише читання)

```powershell
# Python: розрядність має бути 64, шлях — не WindowsApps (Store)
python -c "import sys, struct; print(struct.calcsize('P')*8, sys.executable)"

# Наявні пакети
python -m pip list | Select-String "^(pip|numpy|cffi|bladerf)\s"

# curl / tar (є у Win10 1803+)
(Get-Command curl.exe, tar.exe -ErrorAction SilentlyContinue).Source

# Сторонні libusb (потенційний конфлікт → error 0x7e)
where.exe libusb-1.0.dll

# Версія Windows
(Get-ItemProperty "HKLM:\SOFTWARE\Microsoft\Windows NT\CurrentVersion").DisplayVersion

# Плата: чи бачить Windows (Error до встановлення драйвера — норма)
Get-CimInstance Win32_PnPEntity | Where-Object { $_.DeviceID -match "VID_2CF0|VID_1D50|VID_04B4&PID_00F3" } | Select-Object Status, Name, DeviceID
```

> `Get-PnpDevice` може зависнути — використовуй `Get-CimInstance` (вище).

---

## 2. Python (якщо немає)

| Варіант | Команда / джерело |
| --- | --- |
| winget | `winget install -e --id Python.Python.3.13` (або `3.14`) |
| вручну | [https://www.python.org/downloads/windows/](https://www.python.org/downloads/windows/) → Windows installer (64-bit), ✔ Add to PATH |

Не з Microsoft Store. Після встановлення — новий термінал.

---

## 3. ПЗ bladeRF (драйвер, DLL, CLI, образи FX3/FPGA)

```powershell
cd $env:TEMP
curl.exe -L -o bladeRF-win-installer-2025.10.exe https://github.com/Nuand/bladeRF/releases/download/2025.10/bladeRF-win-installer-2025.10.exe
.\bladeRF-win-installer-2025.10.exe
```

Альтернатива: браузером з [https://github.com/Nuand/bladeRF/releases/tag/2025.10](https://github.com/Nuand/bladeRF/releases/tag/2025.10) → Assets.

- Інсталятор пропонує оновити прошивку плати → плата має бути підключена на цьому етапі.
- Після прошивки помилка `firmware is currently v2.4.0 … requires v2.6.0` — **норма**: потрібен power cycle (від'єднати / під'єднати).
- Образи лишаються на ПК у теці встановлення: `hosted*.rbf`, `fx3_images\bladeRF_fw_v2.6.0.img`. Тека: `Split-Path (Get-Command bladeRF-cli).Source` (образи — на рівень вище від `x64`).

**Перевірка:**

```powershell
bladeRF-cli -p
bladeRF-cli -e version
```

Очікувано: Firmware 2.6.0, FPGA 0.16.0 `(configured by USB host)` — FPGA бібліотека вантажить сама.

Ручна прошивка (якщо треба): `bladeRF-cli -f "<тека>\fx3_images\bladeRF_fw_v2.6.0.img"` → power cycle.

---

## 4. Python-біндинги

`cffi` ставиться автоматично як залежність. `numpy` — окремо, якщо потрібен.

```powershell
cd $env:TEMP
curl.exe -L -o bladeRF-2025.10.zip https://github.com/Nuand/bladeRF/archive/refs/tags/2025.10.zip
tar -xf bladeRF-2025.10.zip
python -m pip install .\bladeRF-2025.10\host\libraries\libbladeRF_bindings\python
```

Варіанти встановлення пакетів:

- **глобально** — `python -m pip install …` (достатньо для пробного запуску);
- **venv** (опційно) — `python -m venv .venv`, далі `.venv\Scripts\python -m pip install …`.

**Перевірка:**

```powershell
bladerf-tool info
python -c "import bladerf; print(bladerf.BladeRF())"
```

---

## 5. Обгортка запуску скрипта (ВІНДОВС)

`D:\Programming\blade\bin\blade_tst_entry.bat`:

```bat
@echo off
setlocal
set "ROOT=%~dp0.."
set "PYTHONPATH=%ROOT%\source;%PYTHONPATH%"
python "%ROOT%\source\blade_tst_entry.py" %*
exit /b %ERRORLEVEL%
```

Виклик: `.\bin\blade_tst_entry.bat --arg1 val "з пробілами"`

---

## 6. Проблеми

| Симптом | Що робити |
| --- | --- |
| `cannot load library 'bladerf.dll': error 0x7e` | Біндинги шукають `bladerf.dll` за іменем → тека `x64` має бути в PATH: `where.exe bladerf.dll`; тимчасово `$env:PATH = "<тека>\x64;$env:PATH"` |
| `python` відкриває Store | Settings → Apps → App execution aliases → вимкнути `python.exe`, `python3.exe` |
| Плата `Status Error` у PnP | Немає драйвера → встановити/перевстановити ПЗ bladeRF |
| `legacy message size` / `firmware update required` | Прошивка \< 2.6.0 або не було power cycle |
| Свій скрипт названо `bladerf.py` | Перейменувати — конфлікт з модулем |