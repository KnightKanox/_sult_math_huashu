# LaTeX 论文撰写 Spec

## Why
当前工作区已经完成问题 1-4 的建模、仿真、验证与结果整理，但还没有把这些内容系统组织成可提交的 LaTeX 论文正文。用户希望基于现有成果撰写论文，且明确要求**不包括第一页**，同时**不包括 AI 相关部分**。

## What Changes
- 生成一份竞赛论文的 **LaTeX 正文稿**，范围从正文开始，不包含第一页内容
- 正文中不单独撰写 AI 工具使用说明、AI 章节或 AI 附录
- 将现有 Q1-Q4 的建模思路、理论推导、仿真结果、图表与结论整合为论文叙述
- 以现有工作区文档与结果文件为唯一事实来源，避免新增未经验证的结论
- 在写作时参考用户提供的论文模板，并将 AI 工具使用章程仅作为排除约束而非正文内容来源

## Impact
- Affected specs: 论文写作、结果汇总、图表组织、结论表述
- Affected code: `version_log/1.0.0/`、`version_log/2.0.0/` 下说明文档与结果文件，`现状描述.md`，`修改日志.md`

## ADDED Requirements
### Requirement: 生成不含第一页的 LaTeX 论文正文
The system SHALL provide a LaTeX manuscript body for the contest paper, excluding the first page content.

#### Scenario: 正文范围明确
- **WHEN** 用户要求撰写 LaTeX 论文且明确说明“不包括第一页”
- **THEN** 生成内容应从正文部分开始，不写第一页封面/标题页/承诺页内容

### Requirement: 正文不包含 AI 相关部分
The system SHALL exclude AI-related sections from the manuscript body.

#### Scenario: 正文排除 AI 内容
- **WHEN** 用户要求论文“不包括 AI 部分”
- **THEN** 生成内容中不得包含 AI 工具使用说明、AI 使用声明、AI 独立章节、AI 附录或围绕 AI 的专门讨论

### Requirement: 论文内容必须基于现有已验证结果
The system SHALL use only validated workspace materials as the factual basis for the manuscript.

#### Scenario: 引用结果与结论
- **WHEN** 论文中需要引用模型结论、参数、图表或数值结果
- **THEN** 内容必须来自现有代码输出、说明文档、版本说明、现状描述或修改日志中的已确认结果

### Requirement: 论文叙述需体现统一建模主线
The system SHALL organize the paper around the established modeling pipeline.

#### Scenario: 组织论文方法部分
- **WHEN** 撰写模型建立与求解方法
- **THEN** 应围绕“几何导通判据 -> 概率云 -> 等效连接体积 -> 渗流解释 -> Monte Carlo 验证 -> 优化求解”展开

### Requirement: 论文需覆盖问题 1-4 的完整成果
The system SHALL include the current project outcomes for all completed problem parts.

#### Scenario: 覆盖题目各问
- **WHEN** 撰写论文正文
- **THEN** 应包含问题 1 的导通判定基线、问题 2 的概率云验证、问题 3 的最低填充率、问题 4 的最低成本组合

### Requirement: 论文措辞需与当前证据强度一致
The system SHALL distinguish between numerical validation and strict proof.

#### Scenario: 表述概率云理论
- **WHEN** 论文讨论概率云理论是否成立
- **THEN** 应表述为“得到充分数值验证/与仿真结果一致”，而不是声称已完成严格数学证明

### Requirement: 论文格式需面向 LaTeX 交付
The system SHALL structure the manuscript in a form directly suitable for LaTeX authoring.

#### Scenario: 输出文稿结构
- **WHEN** 生成论文草稿
- **THEN** 内容应具备章节结构、公式位置、图表引用位置、表格结构与参考文献占位，便于直接落入 LaTeX

## MODIFIED Requirements
### Requirement: 结果文档的使用方式
现有 Markdown 说明文档不再仅作为内部记录，还应作为后续论文正文撰写的直接素材来源。

## REMOVED Requirements
### Requirement: 生成第一页内容
**Reason**: 用户明确要求论文撰写范围不包括第一页。  
**Migration**: 后续仅生成正文 LaTeX 内容，第一页由用户另行处理或后续单独补写。

### Requirement: 生成 AI 相关部分
**Reason**: 用户明确要求正文不包括 AI 部分。  
**Migration**: 后续正文撰写时完全跳过 AI 相关章节、附录和说明。
