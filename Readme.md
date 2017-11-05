主线任务：
基于C++实现一个简单的固定管线软件渲染器

* 渲染一个几何体

* 实现基本的变换（旋转、平移、缩放）

* 实现至少一种光照模型（Blinn-Phong、Lambort等）

* 完成熬测中的Graphics部分

支线任务:

1. 这个任务需要你在你的软件渲染器中实现多线程渲染. 你应该提供一个函数用于设置渲染线程数, 渲染时间应显著少于单线程渲染(存在一个合理的线程数设置, 使得多线程运行时间不超过单线程运行时间的二分之一).
如果你打算做这个任务, 请务必在动手写代码之前就考虑到这一点.

2. 如果你不喜欢多线程, 你可以试试 SIMD (Single Instruction, Multiple Data, 单指令多数据), 尽管 DK 并没有接触过这玩意. 它和GPU的运作方式比较类似, 不过GPU通常在上层做了封装.

3. 在你的软件渲染器中实现光线追踪. 如果你觉得这玩意比较难, 也可以先写一个仅支持镜面的光线追踪.