style:
	@echo "Running code style"
	@pre-commit run --all-files --show-diff-on-failure

install:
	@pip install -r requirements/requirements_build.txt
	@python -m build
	@pip install -q --force-reinstall dist/*.whl

version-info:
	@bash -c "date -u +'Build date: %B %d, %Y %H:%M UTC ShaID: <id>' | xargs -I date sed -i 's/_VERSION_INFO = .*/_VERSION_INFO = \"date\"/g' src/ansys/fluent/core/__init__.py"
	@bash -c "git --no-pager log -n 1 --format='%h' | xargs -I hash sed -i 's/<id>/hash/g' src/ansys/fluent/core/__init__.py"

docker-pull:
	@bash .ci/pull_fluent_image.sh

test-import:
	@python -c "import ansys.fluent.core as pyfluent"

PYTESTEXTRA = --cov=ansys.fluent --cov-report=xml:cov_xml.xml --cov-report=html

unittest: unittest-dev-241

unittest-dev-222:
	@echo "Running unittests"
	@sudo rm -rf /home/ansys/.local/share/ansys_fluent_core/examples
	@pip install -r requirements/requirements_tests.txt
	@python -m pytest --fluent-version=22.2 $(PYTESTEXTRA)

unittest-dev-231:
	@echo "Running unittests"
	@sudo rm -rf /home/ansys/.local/share/ansys_fluent_core/examples
	@pip install -r requirements/requirements_tests.txt
	@python -m pytest --fluent-version=23.1 $(PYTESTEXTRA)

unittest-dev-232:
	@echo "Running unittests"
	@sudo rm -rf /home/ansys/.local/share/ansys_fluent_core/examples
	@pip install -r requirements/requirements_tests.txt
	@python -m pytest --fluent-version=23.2 $(PYTESTEXTRA)

unittest-dev-241:
	@echo "Running unittests"
	@sudo rm -rf /home/ansys/.local/share/ansys_fluent_core/examples
	@pip install -r requirements/requirements_tests.txt
	@python -m pytest --fluent-version=24.1 $(PYTESTEXTRA)

unittest-dev-242:
	@echo "Running unittests"
	@sudo rm -rf /home/ansys/.local/share/ansys_fluent_core/examples
	@pip install -r requirements/requirements_tests.txt
	@python -m pytest --fluent-version=24.2 $(PYTESTEXTRA)

unittest-all-222:
	@echo "Running all unittests"
	@sudo rm -rf /home/ansys/.local/share/ansys_fluent_core/examples
	@pip install -r requirements/requirements_tests.txt
	@python -m pytest --nightly --fluent-version=22.2 $(PYTESTEXTRA)

unittest-all-222-no-codegen:
	@echo "Running all unittests"
	@sudo rm -rf /home/ansys/.local/share/ansys_fluent_core/examples
	@pip install -r requirements/requirements_tests.txt
	@python -m pytest --nightly --fluent-version=22.2 -m "not codegen_required" $(PYTESTEXTRA)

unittest-all-231:
	@echo "Running all unittests"
	@sudo rm -rf /home/ansys/.local/share/ansys_fluent_core/examples
	@pip install -r requirements/requirements_tests.txt
	@python -m pytest --nightly --fluent-version=23.1 $(PYTESTEXTRA)

unittest-all-231-no-codegen:
	@echo "Running all unittests"
	@sudo rm -rf /home/ansys/.local/share/ansys_fluent_core/examples
	@pip install -r requirements/requirements_tests.txt
	@python -m pytest --nightly --fluent-version=23.1 -m "not codegen_required" $(PYTESTEXTRA)

unittest-all-232:
	@echo "Running all unittests"
	@sudo rm -rf /home/ansys/.local/share/ansys_fluent_core/examples
	@pip install -r requirements/requirements_tests.txt
	@python -m pytest --nightly --fluent-version=23.2 $(PYTESTEXTRA)

unittest-all-232-no-codegen:
	@echo "Running all unittests"
	@sudo rm -rf /home/ansys/.local/share/ansys_fluent_core/examples
	@pip install -r requirements/requirements_tests.txt
	@python -m pytest --nightly --fluent-version=23.2 -m "not codegen_required" $(PYTESTEXTRA)

unittest-all-241:
	@echo "Running all unittests"
	@sudo rm -rf /home/ansys/.local/share/ansys_fluent_core/examples
	@pip install -r requirements/requirements_tests.txt
	@python -m pytest --nightly --fluent-version=24.1 $(PYTESTEXTRA)

unittest-all-241-no-codegen:
	@echo "Running all unittests"
	@sudo rm -rf /home/ansys/.local/share/ansys_fluent_core/examples
	@pip install -r requirements/requirements_tests.txt
	@python -m pytest --nightly --fluent-version=24.1 -m "not codegen_required" $(PYTESTEXTRA)

unittest-all-242:
	@echo "Running all unittests"
	@sudo rm -rf /home/ansys/.local/share/ansys_fluent_core/examples
	@pip install -r requirements/requirements_tests.txt
	@python -m pytest --nightly --fluent-version=24.2 $(PYTESTEXTRA)

unittest-solvermode-242:
	@echo "Running all unittests"
	@sudo rm -rf /home/ansys/.local/share/ansys_fluent_core/examples
	@pip install -r requirements/requirements_tests.txt
	@python -m pytest --fluent-version=24.2 --solvermode $(PYTESTEXTRA)

unittest-all-242-no-codegen:
	@echo "Running all unittests"
	@sudo rm -rf /home/ansys/.local/share/ansys_fluent_core/examples
	@pip install -r requirements/requirements_tests.txt
	@python -m pytest --nightly --fluent-version=24.2 -m "not codegen_required" $(PYTESTEXTRA)

api-codegen:
	@echo "Running API codegen"
	@python -m venv env
	@. env/bin/activate
	@pip install -q -e .
	@python codegen/allapigen.py
	@rm -rf env

build-doc-source:
	@sudo rm -rf doc/source/api/meshing/datamodel
	@sudo rm -rf doc/source/api/meshing/tui
	@sudo rm -rf doc/source/api/solver/datamodel
	@sudo rm -rf doc/source/api/solver/tui
	@sudo rm -rf doc/source/api/solver/_autosummary/settings
	@sudo rm -rf /home/ansys/.local/share/ansys_fluent_core/examples
	@pip install -r requirements/requirements_doc.txt
	@xvfb-run make -C doc html

build-all-docs:
	@python doc/settings_rstgen.py
	@sudo rm -rf /home/ansys/.local/share/ansys_fluent_core/examples
	@pip install -r requirements/requirements_doc.txt
	@xvfb-run make -C doc html

compare-flobject:
	@python .ci/compare_flobject.py
