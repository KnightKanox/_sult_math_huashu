# 全局几何与仿真配置常量、边界模式枚举
from enum import Enum

# 微构体盒体半边长（nm）
BOX_HALF = 5000.0
# 微构体盒体边长（nm）
BOX_SIZE = 2.0 * BOX_HALF
# 导通最短距离阈值 δ（nm）
CONNECT_DELTA = 1.8
# 介质A：直圆柱高度（nm）
A_LENGTH = 5000.0
# 介质A：底面半径（nm）
A_RADIUS = 30.0
# 介质B：球体半径（nm）
B_RADIUS = 200.0
# 左电极平面（垂直于X轴，x=-5000）
LEFT_PLANE_X = -BOX_HALF
# 右电极平面（垂直于X轴，x=5000）
RIGHT_PLANE_X = BOX_HALF
# A-A 导通时两轴线段最短距离阈值 = 2*R_A + δ（nm）
AA_SEG_THRESHOLD = 2.0 * A_RADIUS + CONNECT_DELTA
# 介质A(胶囊近似)到电极导通时，轴到平面的距离阈值 = R_A + δ（nm）
AE_AXIS_THRESHOLD = A_RADIUS + CONNECT_DELTA
# 周期重合端点判定容差（nm）
PERIODIC_MATCH_TOL = 1e-6


# 边界模式：决定附件中"行"是否按边界截断规则合并为同一电学导体
class BoundaryMode(Enum):
    # 按题目边界截断规则：越界回绕后仍为同一电学连续导体，跨边界重合片段合并为同一节点
    PERIODIC_CONNECTED = "periodic_connected"
    # 仅几何处理：附件中每一行视为一个独立导体节点，不合并跨边界片段
    WRAPPED_GEOMETRY_ONLY = "wrapped_geometry_only"
