# 🐳 Enhancement Engine Docker Setup

Dieses Dokument beschreibt, wie du das Enhancement Engine Projekt mit Docker und Docker Compose startest.

## 🚀 Quick Start

```bash
# Klone das Repository (falls noch nicht geschehen)
git clone <repository-url>
cd EnhancementEngine

# Starte alle Services
docker compose up -d

# Öffne im Browser
open http://localhost:5000
```

Das war's! 🎉

## 📋 Voraussetzungen

- **Docker**: Version 20.10+
- **Docker Compose**: Version 2.0+

### Installation prüfen
```bash
docker --version
docker compose version
```

## 🏗️ Services

Das Docker Compose Setup beinhaltet folgende Services:

### 1. **enhancement-engine** (Haupt-App)
- **Port**: 5000
- **Beschreibung**: Flask Web-Anwendung
- **URL**: http://localhost:5000

### 2. **redis** (Optional)
- **Port**: 6379  
- **Beschreibung**: Cache für bessere Performance
- **Status**: Automatisch gestartet

### 3. **nginx** (Optional, Production)
- **Port**: 80, 443
- **Beschreibung**: Reverse Proxy für Production
- **Aktivierung**: `docker compose --profile production up -d`

## 🛠️ Kommandos

### Mit Docker Compose
```bash
# Services starten
docker compose up -d

# Services stoppen
docker compose down

# Logs anzeigen
docker compose logs -f

# Einzelner Service
docker compose logs -f enhancement-engine

# Service neustarten
docker compose restart enhancement-engine

# In Container einsteigen
docker compose exec enhancement-engine bash
```

### Mit Makefile (empfohlen)
```bash
# Alle verfügbaren Kommandos anzeigen
make help

# Schnellstart
make up

# Logs anzeigen
make logs

# Services stoppen
make down

# Komplett bereinigen
make clean

# Tests ausführen
make test

# Development Mode
make dev

# Production Mode mit Nginx
make prod
```

## ⚙️ Konfiguration

### Environment Variables

Kopiere `.env.example` zu `.env` und passe die Werte an:

```bash
cp .env.example .env
```

Wichtige Einstellungen:

```env
# Basis-Konfiguration
FLASK_ENV=production
DEMO_EMAIL=demo@example.com
CACHE_ENABLED=true
LOG_LEVEL=INFO

# Sicherheit (für Production)
SECRET_KEY=your-secret-key-here

# NCBI Zugang (optional)
NCBI_EMAIL=your.email@domain.com
NCBI_API_KEY=your-api-key
```

### Volumes

Persistente Daten werden in Docker Volumes gespeichert:

- **enhancement_cache**: Cache-Daten
- **enhancement_sequences**: Sequenz-Cache  
- **enhancement_logs**: Log-Dateien
- **redis_data**: Redis-Daten

## 🔧 Development

### Hot Reload Development
```bash
# Development Mode mit Auto-Reload
make dev

# Oder direkt:
FLASK_ENV=development FLASK_DEBUG=true docker compose up
```

### Code-Änderungen

Template-Dateien werden automatisch gemountet. Für andere Code-Änderungen:

```bash
# Image neu bauen
docker compose build

# Services neustarten  
docker compose up -d
```

### Tests ausführen
```bash
# Alle Tests
make test

# Oder direkt:
docker compose run --rm enhancement-engine python -m pytest tests/ -v

# Linting
make lint
```

## 🌐 Production Deployment

### Mit Nginx Reverse Proxy
```bash
# Production Mode starten
make prod

# Oder direkt:
docker compose --profile production up -d
```

### SSL/HTTPS Setup

1. SSL-Zertifikate in `./ssl/` ablegen:
   ```
   ssl/
   ├── cert.pem
   └── key.pem
   ```

2. Nginx-Konfiguration anpassen (nginx.conf)

3. Mit HTTPS starten:
   ```bash
   docker compose --profile production up -d
   ```

### Health Checks

```bash
# Service-Status prüfen
make health

# Oder direkt:
curl http://localhost:5000/health
```

## 📊 Monitoring

### Logs
```bash
# Alle Services
make logs

# Live Resource-Nutzung
make stats

# Spezifischer Service
docker compose logs -f enhancement-engine
```

### Performance

Die Anwendung ist optimiert für:
- **Caching**: Redis + lokaler Cache
- **Parallelität**: Multi-threaded Flask
- **Resource-Limits**: Definiert in docker-compose.yml

## 🔄 Backup & Restore

### Backup erstellen
```bash
# Automatisches Backup aller Volumes
make backup

# Manuell
docker run --rm -v enhancement_cache:/data -v $(pwd)/backups:/backup alpine tar czf /backup/cache-backup.tar.gz -C /data .
```

### Restore
```bash
# Volume wiederherstellen
docker run --rm -v enhancement_cache:/data -v $(pwd)/backups:/backup alpine tar xzf /backup/cache-backup.tar.gz -C /data
```

## 🐛 Troubleshooting

### Häufige Probleme

**1. Port bereits belegt**
```bash
# Anderen Service auf Port 5000 finden
lsof -i :5000

# Oder anderen Port verwenden
FLASK_PORT=5001 docker compose up -d
```

**2. Permission Errors**
```bash
# Docker-Berechtigungen prüfen
sudo usermod -aG docker $USER
# Neu anmelden erforderlich
```

**3. Out of Memory**
```bash
# Container-Limits erhöhen in docker-compose.yml
services:
  enhancement-engine:
    deploy:
      resources:
        limits:
          memory: 2G
```

**4. Langsame Performance**
```bash
# Cache-Volumes prüfen
docker volume ls
docker volume inspect enhancement_cache

# Redis-Status prüfen
docker compose exec redis redis-cli ping
```

### Debug Mode

```bash
# Verbose Logging
LOG_LEVEL=DEBUG docker compose up

# Container-Zustand prüfen
docker compose ps
docker compose top

# In Container einsteigen
make shell
```

### Komplett-Reset

```bash
# Alles zurücksetzen
make clean

# Images auch löschen
docker compose down --rmi all -v --remove-orphans
docker system prune -a -f
```

## 📱 API Usage

### Web Interface
- **URL**: http://localhost:5000
- **Funktionen**: Gene-Analyse via Web-Form

### REST API
```bash
# Gene analysieren
curl -X POST http://localhost:5000/analyze \
  -H "Content-Type: application/json" \
  -d '{"gene": "COMT", "variant": "Val158Met"}'
```

## 🔐 Sicherheit

### Production Checklist

- [ ] `.env` Datei mit sicheren Werten
- [ ] `SECRET_KEY` gesetzt
- [ ] SSL-Zertifikate konfiguriert  
- [ ] Firewall-Regeln aktiviert
- [ ] Updates regelmäßig durchführen
- [ ] Monitoring eingerichtet
- [ ] Backup-Strategie implementiert

### Security Headers

Nginx ist vorkonfiguriert mit:
- Rate Limiting
- GZIP Compression  
- Security Headers
- SSL/TLS Terminierung

## 🆘 Support

Bei Problemen:

1. **Logs prüfen**: `make logs`
2. **Health Check**: `make health`  
3. **Clean Restart**: `make clean && make up`
4. **GitHub Issues**: Erstelle ein Issue mit Logs

## 📚 Weiterführende Links

- [Docker Compose Dokumentation](https://docs.docker.com/compose/)
- [Flask Deployment Guide](https://flask.palletsprojects.com/en/2.3.x/deploying/)
- [Nginx Configuration](https://nginx.org/en/docs/)
- [Redis Dokumentation](https://redis.io/documentation)