# Render docs/*.svg to high-resolution PNG with headless Chrome/Edge.
#
#   pwsh -File tools\render_diagrams.ps1
#   powershell -NoProfile -ExecutionPolicy Bypass -File tools\render_diagrams.ps1
#
# Output: docs\architecture-code.png, docs\architecture-business.png
#         at 2x device scale (3840 px wide).
#
# NOTE: Chrome writes its progress line to stderr; keep ErrorActionPreference at
# 'Continue' so that harmless native-command noise does not abort the script.
$ErrorActionPreference = 'Continue'

$candidates = @(
    "$env:ProgramFiles\Google\Chrome\Application\chrome.exe",
    "${env:ProgramFiles(x86)}\Google\Chrome\Application\chrome.exe",
    "$env:ProgramFiles\Microsoft\Edge\Application\msedge.exe",
    "${env:ProgramFiles(x86)}\Microsoft\Edge\Application\msedge.exe"
)
$browser = $candidates | Where-Object { Test-Path $_ } | Select-Object -First 1
if (-not $browser) { throw "neither Chrome nor Edge found; install one or render the SVGs by hand" }

$docs = Join-Path (Split-Path -Parent $PSScriptRoot) 'docs'
$scale = 2

foreach ($name in 'architecture-code', 'architecture-business') {
    $svg = Join-Path $docs "$name.svg"
    $png = Join-Path $docs "$name.png"
    if (-not (Test-Path $svg)) { Write-Warning "missing $svg"; continue }

    $head = (Get-Content $svg -TotalCount 3) -join ' '
    $w = [regex]::Match($head, 'width="(\d+)"').Groups[1].Value
    $h = [regex]::Match($head, 'height="(\d+)"').Groups[1].Value
    if (-not $w -or -not $h) { Write-Warning "cannot read size of $svg"; continue }

    Remove-Item $png -ErrorAction SilentlyContinue
    $uri = 'file:///' + ($svg -replace '\\', '/')
    & $browser --headless=new --disable-gpu --hide-scrollbars --no-sandbox `
        --log-level=3 --force-device-scale-factor=$scale --window-size="$w,$h" `
        --default-background-color=FFFFFFFF --screenshot="$png" $uri 2>$null | Out-Null
    if (Test-Path $png) {
        Write-Output ("{0}: {1}x{2} -> {3} ({4:N0} bytes)" -f $name, $w, $h, $png, (Get-Item $png).Length)
    } else {
        Write-Warning "failed to render $name"
    }
}
