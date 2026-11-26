# ChemWand 安装说明

**用话告诉 Claude 画什么，得到能在 ChemDraw 里直接编辑的结构图。**

[English README →](README.md)

ChemWand 是 Claude 的一个扩展。你用中文或英文描述要画的东西——一个分子、一步反应、一条多步路线、
一张底物范围表、一个催化循环——它画好后存成 ChemDraw 的 `.cdxml` 文件，并自动用 ChemDraw 打开，
你可以接着手动修改。每个结构在保存前都会核对一遍，确认画出来的就是你要的那个分子。

全部在你自己的电脑上运行：不用注册、不用 API Key、不用服务器。

---

## 一、准备

只需要两样：

1. **Claude 桌面版**：从 https://claude.com/download 下载，登录你的账号。已经装了的，更新到最新版。
2. **ChemDraw**。

**不需要装 Python 或其他任何软件**，Claude 会自己装好。

## 二、安装（三步，约 2 分钟）

1. 下载安装文件：**[chemwand.mcpb](https://github.com/ClickChemAI/chemwand/releases/latest/download/chemwand.mcpb)**
2. **双击**下载好的 `chemwand.mcpb`，Claude 桌面版会打开一个安装页面。
3. 点 **Install**（安装）。

第一次安装时 Claude 要下载 Python 和化学计算库，需要等一两分钟。装好后就能用了。

> **Install 按钮是灰的、点不了？**
> 说明你登录的是单位（公司、学校）给的 Claude 账号，管理员只允许装批准过的扩展。
> 两个办法：请单位 IT 把 ChemWand 加进允许列表；或者用下面"四、装不了扩展时"的办法。

## 三、开始画图

打开一个普通对话，直接说，例如：

- 用 ChemWand 画一个布洛芬
- 画 Suzuki 偶联：4-溴苯甲醚 + 苯硼酸，Pd(PPh3)4，K2CO3，甲苯/水，90 °C，得到 4-甲氧基联苯，产率 92%
- 做一个底物范围表，R = Me、Et、iPr、Ph，产率分别是 85%、80%、72%、90%
- 画 Buchwald–Hartwig 胺化的催化循环
- 这一步按 0.5 mmol 规模，每种试剂要称多少？

画好的图会自动用 ChemDraw 打开，文件存在个人文件夹里的 `chemwand_figures`：

- Mac：`~/chemwand_figures`
- Windows：`C:\Users\你的用户名\chemwand_figures`

**小提示：** 改了图要重画时，让 Claude 换一个新文件名。同名文件重画，ChemDraw 可能还显示着旧窗口。

## 四、装不了扩展时（单位账号）

1. 在本页点绿色的 **Code** 按钮 → **Download ZIP**，解压，把文件夹放到一个固定位置（之后别再移动）。
2. 打开 Claude 桌面版，切到 **Code** 标签，选择刚才解压的文件夹。
3. 粘贴下面这句发送，中途请求运行命令时点"允许"：

   > Install ChemWand from this folder as a local MCP server, following the README. If Python 3.11 or newer is missing, install it in my user folder with uv. Explain each step to me in Chinese.

4. 装好后**完全退出** Claude（Mac 按 Cmd+Q；Windows 从右下角托盘图标退出），再重新打开。
5. 之后在 **Code** 标签里画图。

## 五、常见问题

| 情况 | 怎么办 |
|---|---|
| Claude 说没有 ChemWand 的工具 | 完全退出 Claude 再打开；还不行，到设置 → 扩展（Extensions）里看 ChemWand 是否已开启 |
| 安装时报错 | 截图报错，到 [Issues](https://github.com/ClickChemAI/chemwand/issues) 提问 |
| 基团没有显示成缩写 | 目前支持 41 种常见缩写（Boc、Cbz、Bn、OTBS、OMe 等），不在表里的会画出完整结构 |
| Windows 上有问题 | Windows 测得比较少，欢迎到 [Issues](https://github.com/ClickChemAI/chemwand/issues) 反馈 |

---

MIT 许可证，可以自由使用、修改、分享。ChemDraw 是 Revvity 公司的商标，ChemWand 是独立项目，与 Revvity 无关，
不包含也不修改任何 ChemDraw 软件。
