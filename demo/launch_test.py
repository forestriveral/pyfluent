import ansys.fluent.core as pyfluent

solver_session = pyfluent.launch_fluent(
    mode="solver",
    version='3d',
    precision='double',
    processor_count=6,
    show_gui=True
    )
solver_session.health_check_service.is_serving