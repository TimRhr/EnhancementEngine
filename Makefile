# Enhancement Engine Docker Management Makefile

.PHONY: help build up down restart logs shell clean test lint dev prod

# Default target
help:
	@echo "Enhancement Engine Docker Commands:"
	@echo "  build     - Build Docker images"
	@echo "  up        - Start all services"
	@echo "  down      - Stop all services"
	@echo "  restart   - Restart all services"
	@echo "  logs      - Show logs"
	@echo "  shell     - Get shell in main container"
	@echo "  clean     - Clean up containers and volumes"
	@echo "  test      - Run tests in container"
	@echo "  lint      - Run linting in container"
	@echo "  dev       - Start in development mode"
	@echo "  prod      - Start in production mode with Nginx"
	@echo ""
	@echo "Quick start: make up"

# Build Docker images
build:
	@echo "Building Enhancement Engine Docker images..."
	docker compose build

# Start all services
up:
	@echo "Starting Enhancement Engine..."
	docker compose up -d
	@echo "Enhancement Engine is running at http://localhost:5000"

# Stop all services
down:
	@echo "Stopping Enhancement Engine..."
	docker compose down

# Restart all services
restart: down up

# Show logs
logs:
	docker compose logs -f

# Get shell in main container
shell:
	docker compose exec enhancement-engine bash

# Clean up everything
clean:
	@echo "Cleaning up Enhancement Engine containers and volumes..."
	docker compose down -v --remove-orphans
	docker system prune -f
	@echo "Cleanup complete!"

# Run tests
test:
	@echo "Running tests..."
	docker compose run --rm enhancement-engine python -m pytest tests/ -v

# Run linting
lint:
	@echo "Running linting..."
	docker compose run --rm enhancement-engine python -m flake8 enhancement_engine/
	docker compose run --rm enhancement-engine python -m black --check enhancement_engine/

# Development mode (with hot reload)
dev:
	@echo "Starting in development mode..."
	FLASK_ENV=development FLASK_DEBUG=true docker compose up

# Production mode (with Nginx)
prod:
	@echo "Starting in production mode with Nginx..."
	docker compose --profile production up -d
	@echo "Enhancement Engine is running at http://localhost (Nginx proxy)"

# Initialize environment file
init:
	@if [ ! -f .env ]; then \
		echo "Creating .env file from .env.example..."; \
		cp .env.example .env; \
		echo "Please edit .env file with your configuration"; \
	else \
		echo ".env file already exists"; \
	fi

# Check service health
health:
	@echo "Checking service health..."
	@docker compose ps
	@echo ""
	@curl -f http://localhost:5000/health 2>/dev/null && echo "✅ App is healthy" || echo "❌ App is not responding"

# View real-time resource usage
stats:
	docker stats

# Backup data volumes
backup:
	@echo "Creating backup of data volumes..."
	@mkdir -p backups
	docker run --rm -v enhancement_cache:/data -v $(PWD)/backups:/backup alpine tar czf /backup/cache-$(shell date +%Y%m%d-%H%M%S).tar.gz -C /data .
	docker run --rm -v enhancement_sequences:/data -v $(PWD)/backups:/backup alpine tar czf /backup/sequences-$(shell date +%Y%m%d-%H%M%S).tar.gz -C /data .
	@echo "Backup complete! Files saved in ./backups/"

# Update dependencies
update:
	@echo "Updating dependencies..."
	docker compose build --no-cache
	@echo "Dependencies updated!"