# 随机生成层：介质A圆柱的随机方向/位置采样与周期回绕切段生成
# 本包用于问题 2：随机微构体生成（球面均匀方向、盒内均匀位置、跨边界圆柱解析切段）
from .random_orientation import random_unit_vector
from .random_position import random_position_in_box
from .medium_generator import (generate_a_cylinder, generate_batch,
                               make_cylinder_segments)
