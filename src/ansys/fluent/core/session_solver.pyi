from ansys.fluent.core.datamodel_241.preferences import Root as preferences_root
from ansys.fluent.core.datamodel_241.solverworkflow import Root as solverworkflow_root
from ansys.fluent.core.datamodel_241.workflow import Root as workflow_root
from ansys.fluent.core.solver.settings_241.current_parametric_study import (
    current_parametric_study,
)
from ansys.fluent.core.solver.settings_241.file import file
from ansys.fluent.core.solver.settings_241.mesh import mesh
from ansys.fluent.core.solver.settings_241.parallel import parallel
from ansys.fluent.core.solver.settings_241.parametric_studies import parametric_studies
from ansys.fluent.core.solver.settings_241.report import report
from ansys.fluent.core.solver.settings_241.results import results
from ansys.fluent.core.solver.settings_241.server import server
from ansys.fluent.core.solver.settings_241.setup import setup
from ansys.fluent.core.solver.settings_241.solution import solution
from ansys.fluent.core.solver.tui_241 import main_menu
from ansys.fluent.core.systemcoupling import SystemCoupling

class Solver:
    def build_from_fluent_connection(self, fluent_connection): ...
    @property
    def version(self): ...
    @property
    def tui(self) -> main_menu: ...
    @property
    def workflow(self) -> workflow_root: ...
    @property
    def system_coupling(self) -> SystemCoupling: ...
    @property
    def preferences(self) -> preferences_root: ...
    def read_case_lightweight(self, file_name: str): ...
    def read_case(self, file_name: str): ...
    def write_case(self, file_name: str): ...
    @property
    def file(self) -> file: ...
    @property
    def mesh(self) -> mesh: ...
    @property
    def server(self) -> server: ...
    @property
    def setup(self) -> setup: ...
    @property
    def solution(self) -> solution: ...
    @property
    def results(self) -> results: ...
    @property
    def design(self) -> design: ...
    @property
    def parametric_studies(self) -> parametric_studies: ...
    @property
    def current_parametric_study(self) -> current_parametric_study: ...
    @property
    def parameters(self) -> parameters: ...
    @property
    def parallel(self) -> parallel: ...
