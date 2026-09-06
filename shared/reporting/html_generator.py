"""Standalone HTML report generator using Jinja2 with modern responsive styling."""

from pathlib import Path
from typing import Any, Dict, List, Optional

import jinja2

REPORT_HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>{{ title }}</title>
  <style>
    :root {
      --primary: #1e3a8a;
      --primary-light: #3b82f6;
      --secondary: #0f766e;
      --bg: #f8fafc;
      --surface: #ffffff;
      --text: #0f172a;
      --text-muted: #64748b;
      --border: #e2e8f0;
      --accent: #f59e0b;
    }
    * { box-sizing: border-box; margin: 0; padding: 0; }
    body {
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
      background-color: var(--bg);
      color: var(--text);
      line-height: 1.6;
      padding: 2rem 1rem;
    }
    .container { max-width: 1100px; margin: 0 auto; }
    header {
      background: linear-gradient(135deg, var(--primary), var(--secondary));
      color: white;
      padding: 2.5rem 2rem;
      border-radius: 12px;
      margin-bottom: 2rem;
      box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1);
    }
    header h1 { font-size: 2rem; margin-bottom: 0.5rem; }
    header p { opacity: 0.9; font-size: 1rem; }
    .badge {
      display: inline-block;
      background: rgba(255, 255, 255, 0.2);
      padding: 0.2rem 0.6rem;
      border-radius: 6px;
      font-size: 0.85rem;
      margin-top: 0.5rem;
    }
    .grid-stats {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
      gap: 1rem;
      margin-bottom: 2rem;
    }
    .stat-card {
      background: var(--surface);
      border: 1px solid var(--border);
      border-radius: 10px;
      padding: 1.25rem;
      box-shadow: 0 1px 3px rgba(0,0,0,0.05);
    }
    .stat-card .label { font-size: 0.8rem; color: var(--text-muted); text-transform: uppercase; letter-spacing: 0.5px; }
    .stat-card .value { font-size: 1.5rem; font-weight: 700; color: var(--primary); margin-top: 0.25rem; }
    .card {
      background: var(--surface);
      border: 1px solid var(--border);
      border-radius: 10px;
      padding: 1.5rem;
      margin-bottom: 2rem;
      box-shadow: 0 1px 3px rgba(0,0,0,0.05);
    }
    .card h2 {
      font-size: 1.25rem;
      color: var(--primary);
      margin-bottom: 1rem;
      border-bottom: 2px solid var(--border);
      padding-bottom: 0.5rem;
    }
    table { width: 100%; border-collapse: collapse; margin-top: 1rem; font-size: 0.9rem; }
    th, td { padding: 0.75rem 1rem; text-align: left; border-bottom: 1px solid var(--border); }
    th { background: #f1f5f9; font-weight: 600; color: var(--text-muted); }
    tr:hover { background-color: #f8fafc; }
    .figure-container { text-align: center; margin: 1.5rem 0; }
    .figure-container img { max-width: 100%; height: auto; border-radius: 8px; border: 1px solid var(--border); }
    .figure-caption { font-size: 0.85rem; color: var(--text-muted); margin-top: 0.5rem; }
    .disclaimer {
      background-color: #fffbeb;
      border-left: 4px solid var(--accent);
      padding: 1rem;
      border-radius: 6px;
      font-size: 0.9rem;
      color: #92400e;
      margin-top: 2rem;
    }
    footer { text-align: center; color: var(--text-muted); font-size: 0.85rem; margin-top: 3rem; }
  </style>
</head>
<body>
  <div class="container">
    <header>
      <h1>{{ title }}</h1>
      <p>{{ subtitle }}</p>
      <div class="badge">Generated: {{ timestamp }} | Pipeline: {{ pipeline_name }}</div>
    </header>

    {% if stats %}
    <div class="grid-stats">
      {% for stat in stats %}
      <div class="stat-card">
        <div class="label">{{ stat.label }}</div>
        <div class="value">{{ stat.value }}</div>
      </div>
      {% endfor %}
    </div>
    {% endif %}

    {% for section in sections %}
    <div class="card">
      <h2>{{ section.title }}</h2>
      {% if section.content %}
      <p>{{ section.content | safe }}</p>
      {% endif %}

      {% if section.figure_url %}
      <div class="figure-container">
        <img src="{{ section.figure_url }}" alt="{{ section.title }}">
        {% if section.figure_caption %}
        <div class="figure-caption">{{ section.figure_caption }}</div>
        {% endif %}
      </div>
      {% endif %}

      {% if section.table_data %}
      <div style="overflow-x: auto;">
        <table>
          <thead>
            <tr>
              {% for col in section.table_columns %}
              <th>{{ col }}</th>
              {% endfor %}
            </tr>
          </thead>
          <tbody>
            {% for row in section.table_data %}
            <tr>
              {% for cell in row %}
              <td>{{ cell }}</td>
              {% endfor %}
            </tr>
            {% endfor %}
          </tbody>
        </table>
      </div>
      {% endif %}
    </div>
    {% endfor %}

    {% if disclaimer %}
    <div class="disclaimer">
      <strong>Scientific & Ethical Notice:</strong> {{ disclaimer }}
    </div>
    {% endif %}

    <footer>
      <p>Computational Biology & Bioinformatics Engineering Portfolio • Built with Python 3.11+</p>
    </footer>
  </div>
</body>
</html>
"""


def render_html_report(
    title: str,
    subtitle: str,
    pipeline_name: str,
    timestamp: str,
    stats: Optional[List[Dict[str, str]]] = None,
    sections: Optional[List[Dict[str, Any]]] = None,
    disclaimer: Optional[str] = None,
    output_path: str | Path = "report.html",
) -> Path:
    """Render a structured report dictionary into a clean, standalone HTML file."""
    template = jinja2.Template(REPORT_HTML_TEMPLATE)
    html_out = template.render(
        title=title,
        subtitle=subtitle,
        pipeline_name=pipeline_name,
        timestamp=timestamp,
        stats=stats or [],
        sections=sections or [],
        disclaimer=disclaimer,
    )

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        f.write(html_out)

    return out
